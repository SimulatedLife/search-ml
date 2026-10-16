#!/usr/bin/env python3
"""
searchml-hpc : the HPC-side job manager used by the SEARCH-ML web portal.

The web server talks to the cluster over SSH and only ever runs this program.
Every command prints one JSON object on stdout (except ``get``, which streams a
file) and exits 0 on success, 1 on error.

Commands
--------
  ping                         server health, scheduler, installed libraries
  prepare <job_id>             create the job directory
  put <job_id> <file>          store stdin as <file> in the job directory
  submit <job_id>              queue the job (local / SLURM / PBS)
  status <job_id>              state, stage, progress, ETA, errors, result files
  get <job_id> <file>          stream a result/log file to stdout
  cancel <job_id>              cancel a queued or running job
  list [--limit N]             recent jobs (admin)
  cleanup [--days N] [--dry-run]   delete finished jobs older than N days
  selftest [--simulate|--real] run a test job end to end
  version

Configuration: hpc/searchml-hpc.conf (KEY=VALUE, see searchml-hpc.conf.example),
or the file named by $SEARCHML_HPC_CONF.

Only the Python standard library is used (Python >= 3.6).
"""

import fcntl
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime

VERSION = "2.0.0"
HERE = os.path.dirname(os.path.realpath(__file__))
HPC_DIR = os.path.dirname(HERE)
REPO_DIR = os.path.dirname(HPC_DIR)

JOB_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{5,63}$")
UPLOAD_NAMES = re.compile(
    r"^(job\.json|receptor\.pdb|custom_molecules\.sdf|"
    r"(single_molecule|query_ligand)\.(sdf|mol|mol2|pdb|smi))$")
DOWNLOAD_NAMES = {
    "results.csv", "results.csv.gz", "results_preview.csv", "summary.json", "job.log", "error.log",
    "state.json", "progress.json", "prepare.json", "parameter.txt", "PROT_LIG.pdb", "leap.log",
    "job.json", "Ligand.sdf",
}
TERMINAL = {"COMPLETED", "FAILED", "CANCELLED"}

DEFAULTS = {
    "SEARCH_ML_HOME": REPO_DIR,
    "JOBS_DIR": os.path.join(os.path.expanduser("~"), "searchml-jobs"),
    "DATASETS_DIR": "",
    "MODELS_DIR": "",
    "CONDA_BASE": "",
    "CONDA_ENV": "search-ml",
    "PYTHON": "",
    "SCHEDULER": "local",
    "MAX_LOCAL_JOBS": "2",
    "SMALL_CPUS": "4",
    "SMALL_MEM": "16G",
    "SMALL_TIME": "02:00:00",
    "LARGE_CPUS": "16",
    "LARGE_MEM": "96G",
    "LARGE_TIME": "24:00:00",
    "PARTITION": "",
    "ACCOUNT": "",
    "QOS": "",
    "EXTRA_SBATCH": "",
    "EXTRA_QSUB": "",
    "MODULES": "",
    "RESULT_LIMIT_LARGE": "100000",
    "PREVIEW_ROWS": "2000",
    "MAX_UPLOAD_MB": "50",
    "MAX_CUSTOM_MOLECULES": "5000",
    "MAX_ACTIVE_JOBS": "50",
    "JOB_RETENTION_DAYS": "30",
    "GZIP_RESULTS_MB": "50",
}


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

class HpcError(Exception):
    pass


def out(obj, code=0):
    sys.stdout.write(json.dumps(obj, separators=(",", ":"), default=str) + "\n")
    sys.stdout.flush()
    sys.exit(code)


def fail(msg, **extra):
    payload = {"ok": False, "error": msg}
    payload.update(extra)
    out(payload, 1)


def now_iso():
    return datetime.now().astimezone().isoformat(timespec="seconds")


def read_json(path, default=None):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return default


def write_json(path, payload):
    tmp = "%s.tmp%d" % (path, os.getpid())
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def tail(path, n=20, max_bytes=65536):
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - max_bytes))
            data = fh.read().decode("utf-8", "replace")
        return "\n".join(data.splitlines()[-n:])
    except OSError:
        return ""


def load_config():
    cfg = dict(DEFAULTS)
    path = os.environ.get("SEARCHML_HPC_CONF") or os.path.join(HPC_DIR, "searchml-hpc.conf")
    if os.path.isfile(path):
        with open(path, "r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                key, _, val = line.partition("=")
                key = key.strip()
                val = val.split(" #", 1)[0].strip()
                if len(val) >= 2 and val[0] == val[-1] and val[0] in "\"'":
                    val = val[1:-1]
                cfg[key] = os.path.expandvars(os.path.expanduser(val))
    cfg["_CONF_PATH"] = path if os.path.isfile(path) else ""
    if not cfg["PYTHON"]:
        if cfg["CONDA_BASE"]:
            cand = os.path.join(cfg["CONDA_BASE"], "envs", cfg["CONDA_ENV"], "bin", "python")
            cfg["PYTHON"] = cand if os.path.isfile(cand) else sys.executable
        else:
            cfg["PYTHON"] = sys.executable
    cfg["ENV_PREFIX"] = os.path.dirname(os.path.dirname(os.path.realpath(cfg["PYTHON"])))
    if not cfg["DATASETS_DIR"]:
        cfg["DATASETS_DIR"] = os.path.join(cfg["SEARCH_ML_HOME"], "datasets")
    if not cfg["MODELS_DIR"]:
        cfg["MODELS_DIR"] = os.path.join(cfg["SEARCH_ML_HOME"], "models")
    cfg["SCHEDULER"] = cfg["SCHEDULER"].lower()
    return cfg


def ival(cfg, key):
    try:
        return int(str(cfg.get(key, "0")).strip() or 0)
    except ValueError:
        return int(DEFAULTS.get(key, "0") or 0)


def job_dir(cfg, job_id, must_exist=True):
    if not job_id or not JOB_ID_RE.match(job_id):
        raise HpcError("invalid job id")
    path = os.path.join(cfg["JOBS_DIR"], job_id)
    if must_exist and not os.path.isdir(path):
        raise HpcError("job not found")
    return path


def pipeline_env(cfg, cpus=None):
    env = os.environ.copy()
    prefix = cfg["ENV_PREFIX"]
    env["PATH"] = os.path.join(prefix, "bin") + os.pathsep + env.get("PATH", "")
    env.setdefault("AMBERHOME", prefix)
    if not os.path.isdir(os.path.join(env["AMBERHOME"], "dat", "leap")):
        env["AMBERHOME"] = prefix
    env["CONDA_PREFIX"] = prefix
    env["SEARCH_ML_HOME"] = cfg["SEARCH_ML_HOME"]
    env["SEARCHML_DATASETS"] = cfg["DATASETS_DIR"]
    env["SEARCHML_MODELS"] = cfg["MODELS_DIR"]
    env["SEARCHML_GZIP_MB"] = str(cfg.get("GZIP_RESULTS_MB", "50"))
    for k in ("SEARCHML_DATABASES_TSV", "SEARCHML_LIBRARIES_TSV"):
        if cfg.get(k):
            env[k] = cfg[k]
    if cpus:
        env["SEARCHML_CPUS"] = str(cpus)
    env["PYTHONUNBUFFERED"] = "1"
    return env


def run(cmd, timeout=20):
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        o, e = p.communicate(timeout=timeout)
        return p.returncode, o, e
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 127, "", str(exc)


# --------------------------------------------------------------------------
# catalogs (same files the pipeline reads)
# --------------------------------------------------------------------------

def catalog(cfg, kind):
    env_key = "SEARCHML_DATABASES_TSV" if kind == "databases" else "SEARCHML_LIBRARIES_TSV"
    path = cfg.get(env_key) or os.path.join(cfg["SEARCH_ML_HOME"], "config", kind + ".tsv")
    entries = []
    if not os.path.isfile(path):
        return entries
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            entries.append({"key": cols[0].strip().lower(), "aliases": [a for a in cols[1].split(",") if a],
                            "label": cols[2].strip(), "path": cols[3].strip(),
                            "size": (cols[4].strip() if len(cols) > 4 else "small") or "small"})
    return entries


def dataset_info(cfg, entry, kind):
    path = entry["path"]
    if not os.path.isabs(path):
        path = os.path.join(cfg["DATASETS_DIR"], path)
    main = "dataset.csv" if kind == "databases" else "protein.csv"
    info = {"label": entry["label"], "size": entry["size"], "available": False, "rows": None}
    files = []
    if os.path.isdir(path):
        for cand in (main, main + ".gz"):
            if os.path.isfile(os.path.join(path, cand)):
                files = [os.path.join(path, cand)]
                break
        if not files:
            parts = os.path.join(path, "parts") if os.path.isdir(os.path.join(path, "parts")) else path
            files = [os.path.join(parts, f) for f in os.listdir(parts) if f.endswith((".csv", ".csv.gz"))]
        cnt = os.path.join(path, "count.txt")
        if os.path.isfile(cnt):
            try:
                info["rows"] = int(open(cnt).read().split()[0])
            except (ValueError, IndexError):
                pass
    elif os.path.isfile(path):
        files = [path]
    if files:
        info["available"] = True
        if info["rows"] is None and len(files) == 1 and not files[0].endswith(".gz") \
                and os.path.getsize(files[0]) < 64 * 1024 * 1024:
            with open(files[0], "rb") as fh:
                info["rows"] = max(0, sum(chunk.count(b"\n") for chunk in iter(lambda: fh.read(1 << 20), b"")) - 1)
    return info


# --------------------------------------------------------------------------
# scheduler back-ends
# --------------------------------------------------------------------------

def job_profile(cfg, spec):
    """Return 'large' or 'small' for resource selection."""
    if (spec.get("mode") or "") == "forward":
        f = spec.get("forward") or {}
        if (f.get("source") or "database") == "database":
            key = str(f.get("database", "")).lower()
            for e in catalog(cfg, "databases"):
                if e["key"] == key or key in [a.lower() for a in e["aliases"]]:
                    return e["size"]
    else:
        key = str((spec.get("reverse") or {}).get("library", "homo")).lower()
        for e in catalog(cfg, "libraries"):
            if e["key"] == key:
                return e["size"]
    return "small"


def render_template(name, values):
    with open(os.path.join(HPC_DIR, "templates", name), "r", encoding="utf-8") as fh:
        text = fh.read()
    for k, v in values.items():
        text = text.replace("@@%s@@" % k, str(v))
    return text


def template_values(cfg, jdir, profile):
    p = "LARGE" if profile == "large" else "SMALL"
    modules = " ".join(shlex.quote(m) for m in cfg.get("MODULES", "").split())
    return {
        "JOB_DIR": shlex.quote(jdir),
        "JOB_NAME": "sml-" + os.path.basename(jdir)[-20:],
        "CPUS": ival(cfg, p + "_CPUS") or 1,
        "MEM": cfg.get(p + "_MEM", "16G"),
        "TIME": cfg.get(p + "_TIME", "02:00:00"),
        "PYTHON": shlex.quote(cfg["PYTHON"]),
        "ENV_PREFIX": shlex.quote(cfg["ENV_PREFIX"]),
        "CONDA_BASE": shlex.quote(cfg.get("CONDA_BASE", "")),
        "CONDA_ENV": shlex.quote(cfg.get("CONDA_ENV", "search-ml")),
        "SEARCH_ML_HOME": shlex.quote(cfg["SEARCH_ML_HOME"]),
        "DATASETS_DIR": shlex.quote(cfg["DATASETS_DIR"]),
        "MODELS_DIR": shlex.quote(cfg["MODELS_DIR"]),
        "GZIP_MB": cfg.get("GZIP_RESULTS_MB", "50"),
        "MODULE_LINE": ("module load " + modules) if modules else ": # no environment modules configured",
        "PIPELINE": shlex.quote(os.path.join(cfg["SEARCH_ML_HOME"], "scripts", "pipeline.py")),
    }


def submit_slurm(cfg, jdir, profile):
    vals = template_values(cfg, jdir, profile)
    script = os.path.join(jdir, "submit.slurm")
    with open(script, "w") as fh:
        fh.write(render_template("slurm.sh", vals))
    cmd = ["sbatch", "--parsable", "--job-name", vals["JOB_NAME"], "--chdir", jdir,
           "--output", os.path.join(jdir, "scheduler.out"), "--cpus-per-task", str(vals["CPUS"]),
           "--mem", vals["MEM"], "--time", vals["TIME"], "--nodes", "1", "--ntasks", "1"]
    for key, flag in (("PARTITION", "--partition"), ("ACCOUNT", "--account"), ("QOS", "--qos")):
        if cfg.get(key):
            cmd += [flag, cfg[key]]
    cmd += shlex.split(cfg.get("EXTRA_SBATCH", ""))
    cmd.append(script)
    rc, o, e = run(cmd, timeout=60)
    if rc != 0:
        raise HpcError("sbatch failed: " + (e.strip() or o.strip())[:400])
    return o.strip().split(";")[0]


def submit_pbs(cfg, jdir, profile):
    vals = template_values(cfg, jdir, profile)
    script = os.path.join(jdir, "submit.pbs")
    with open(script, "w") as fh:
        fh.write(render_template("pbs.sh", vals))
    cmd = ["qsub", "-N", vals["JOB_NAME"][:15], "-j", "oe", "-o", os.path.join(jdir, "scheduler.out"),
           "-l", "select=1:ncpus=%s:mem=%s" % (vals["CPUS"], vals["MEM"].lower().replace("g", "gb")
                                               if vals["MEM"][-1:].lower() == "g" else vals["MEM"]),
           "-l", "walltime=%s" % vals["TIME"]]
    if cfg.get("PARTITION"):
        cmd += ["-q", cfg["PARTITION"]]
    if cfg.get("ACCOUNT"):
        cmd += ["-A", cfg["ACCOUNT"]]
    cmd += shlex.split(cfg.get("EXTRA_QSUB", ""))
    cmd.append(script)
    rc, o, e = run(cmd, timeout=60)
    if rc != 0:
        raise HpcError("qsub failed: " + (e.strip() or o.strip())[:400])
    return o.strip().splitlines()[-1] if o.strip() else ""


def submit_local(cfg, jdir, profile):
    log = open(os.path.join(jdir, "scheduler.out"), "ab")
    cmd = [sys.executable, os.path.realpath(__file__), "_local_worker", os.path.basename(jdir), profile]
    env = os.environ.copy()
    if cfg.get("_CONF_PATH"):
        env["SEARCHML_HPC_CONF"] = cfg["_CONF_PATH"]
    p = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=log, stderr=log, cwd=jdir, env=env,
                         start_new_session=True, close_fds=True)
    log.close()
    return "local-%d" % p.pid


def local_worker(cfg, job_id, profile):
    """Runs detached: waits for a free slot, then executes the pipeline."""
    jdir = job_dir(cfg, job_id)
    with open(os.path.join(jdir, "worker.pid"), "w") as fh:
        fh.write("%d\n" % os.getpid())
    slots_dir = os.path.join(cfg["JOBS_DIR"], ".slots")
    os.makedirs(slots_dir, exist_ok=True)
    n_slots = max(1, ival(cfg, "MAX_LOCAL_JOBS"))
    lock_fh = None
    while lock_fh is None:
        if os.path.exists(os.path.join(jdir, "CANCEL_REQUESTED")):
            mark_cancelled(jdir)
            return 0
        for i in range(1, n_slots + 1):
            fh = open(os.path.join(slots_dir, "slot_%d.lock" % i), "w")
            try:
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
                lock_fh = fh
                break
            except OSError:
                fh.close()
        if lock_fh is None:
            time.sleep(5)
    cpus = ival(cfg, "LARGE_CPUS" if profile == "large" else "SMALL_CPUS") or 1
    cmd = [cfg["PYTHON"], os.path.join(cfg["SEARCH_ML_HOME"], "scripts", "pipeline.py"), "run", jdir, "--quiet"]
    child = subprocess.Popen(cmd, cwd=jdir, env=pipeline_env(cfg, cpus))

    def forward(signum, _frame):
        try:
            child.send_signal(signum)
        except OSError:
            pass
    signal.signal(signal.SIGTERM, forward)
    signal.signal(signal.SIGINT, forward)
    rc = child.wait()
    lock_fh.close()
    state = (read_json(os.path.join(jdir, "state.json"), {}) or {}).get("state")
    if state not in TERMINAL:  # pipeline died without recording a state (e.g. killed -9)
        mark_failed(jdir, "The job process exited unexpectedly (code %s)." % rc)
    return rc


def mark_failed(jdir, message):
    st = read_json(os.path.join(jdir, "state.json"), {}) or {}
    st.update(state="FAILED", message=message, finished_at=now_iso(), updated_at=now_iso())
    write_json(os.path.join(jdir, "state.json"), st)
    if not os.path.exists(os.path.join(jdir, "error.log")):
        with open(os.path.join(jdir, "error.log"), "w") as fh:
            fh.write(message + "\n")
            sched = tail(os.path.join(jdir, "scheduler.out"), 30)
            if sched:
                fh.write("\n--- Details ---\n" + sched + "\n")
    for f in ("QUEUED", "RUNNING"):
        try:
            os.remove(os.path.join(jdir, f))
        except FileNotFoundError:
            pass
    open(os.path.join(jdir, "FAILED"), "a").close()


def mark_cancelled(jdir):
    st = read_json(os.path.join(jdir, "state.json"), {}) or {}
    st.update(state="CANCELLED", finished_at=now_iso(), updated_at=now_iso(), message="Cancelled by user")
    write_json(os.path.join(jdir, "state.json"), st)
    for f in ("QUEUED", "RUNNING"):
        try:
            os.remove(os.path.join(jdir, f))
        except FileNotFoundError:
            pass
    open(os.path.join(jdir, "CANCELLED"), "a").close()


def scheduler_state(cfg, sched, sid):
    """Return (state, reason) from the scheduler; state None if the job is unknown/finished."""
    if not sid:
        return None, ""
    if sched == "slurm":
        rc, o, _ = run(["squeue", "-h", "-j", sid, "-o", "%T|%r"], timeout=15)
        if rc == 0 and o.strip():
            parts = o.strip().splitlines()[0].split("|")
            return parts[0], (parts[1] if len(parts) > 1 else "")
        rc, o, _ = run(["sacct", "-j", sid, "-X", "-n", "-P", "-o", "State"], timeout=15)
        if rc == 0 and o.strip():
            return "ENDED:" + o.strip().splitlines()[0].split()[0], ""
        return None, ""
    if sched == "pbs":
        rc, o, _ = run(["qstat", "-f", sid], timeout=15)
        if rc == 0:
            m = re.search(r"job_state\s*=\s*(\w)", o)
            code = m.group(1) if m else "?"
            return {"Q": "PENDING", "R": "RUNNING", "H": "HELD", "E": "EXITING", "F": "ENDED:DONE"}.get(code, code), ""
        return None, ""
    if sched == "local":
        pid = sid.replace("local-", "")
        if pid.isdigit():
            return ("RUNNING", "") if local_pid_alive(int(pid)) else (None, "")
    return None, ""


def local_pid_alive(pid):
    """True if the local worker process exists and is not a zombie / a recycled PID."""
    proc = "/proc/%d" % pid
    if os.path.isdir("/proc/self"):
        try:
            with open(proc + "/stat", "r") as fh:
                state = fh.read().rsplit(")", 1)[1].split()[0]
            if state in ("Z", "X"):
                return False
            with open(proc + "/cmdline", "rb") as fh:
                return b"_local_worker" in fh.read()
        except (OSError, IndexError):
            return False
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


SLURM_END_MESSAGES = {
    "TIMEOUT": "The job reached the cluster time limit before finishing.",
    "OUT_OF_MEMORY": "The job ran out of memory on the compute node.",
    "NODE_FAIL": "The compute node failed while running the job.",
    "CANCELLED": "The job was cancelled on the cluster.",
    "FAILED": "The job failed on the cluster.",
    "PREEMPTED": "The job was pre-empted by the scheduler.",
}


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------

def cmd_ping(cfg, args):
    deep = "--deep" in args
    os.makedirs(cfg["JOBS_DIR"], exist_ok=True)
    st = os.statvfs(cfg["JOBS_DIR"])
    dbs = {e["key"]: dataset_info(cfg, e, "databases") for e in catalog(cfg, "databases")}
    libs = {e["key"]: dataset_info(cfg, e, "libraries") for e in catalog(cfg, "libraries")}
    models = os.path.join(cfg["MODELS_DIR"])
    model_files = ["scaler_final.joblib", "final_meta_model.joblib", "final_xgb_model.joblib",
                   "final_cat_model.joblib", "final_rf_model.joblib", "final_lgbm_model.joblib"]
    missing_models = [m for m in model_files if not os.path.isfile(os.path.join(models, m))]
    tleap = os.path.join(cfg["ENV_PREFIX"], "bin", "tleap")
    active = 0
    for name in os.listdir(cfg["JOBS_DIR"]):
        if os.path.exists(os.path.join(cfg["JOBS_DIR"], name, "RUNNING")) or \
                os.path.exists(os.path.join(cfg["JOBS_DIR"], name, "QUEUED")):
            active += 1
    problems = []
    if missing_models:
        problems.append("missing models: " + ", ".join(missing_models))
    if not os.path.isfile(cfg["PYTHON"]):
        problems.append("python not found: " + cfg["PYTHON"])
    if not os.path.isfile(tleap) and not shutil.which("tleap"):
        problems.append("tleap not found in " + cfg["ENV_PREFIX"])
    if deep:
        rc, o, e = run([cfg["PYTHON"], os.path.join(cfg["SEARCH_ML_HOME"], "scripts", "pipeline.py"), "check"],
                       timeout=180)
        if rc != 0:
            problems.append("pipeline check failed: " + (o + e).strip()[-600:])
    out({
        "ok": True,
        "version": VERSION,
        "host": os.uname()[1],
        "time": now_iso(),
        "scheduler": cfg["SCHEDULER"],
        "jobs_dir": cfg["JOBS_DIR"],
        "disk_free_gb": round(st.f_bavail * st.f_frsize / 1e9, 1),
        "active_jobs": active,
        "accepting": active < max(1, ival(cfg, "MAX_ACTIVE_JOBS")),
        "ready": not problems,
        "problems": problems,
        "databases": dbs,
        "libraries": libs,
        "limits": {"max_upload_mb": ival(cfg, "MAX_UPLOAD_MB"),
                   "max_custom_molecules": ival(cfg, "MAX_CUSTOM_MOLECULES"),
                   "result_limit_large": ival(cfg, "RESULT_LIMIT_LARGE"),
                   "retention_days": ival(cfg, "JOB_RETENTION_DAYS")},
    })


def cmd_prepare(cfg, args):
    if len(args) != 1:
        raise HpcError("usage: prepare <job_id>")
    jdir = job_dir(cfg, args[0], must_exist=False)
    if os.path.exists(os.path.join(jdir, "state.json")):
        raise HpcError("job already exists")
    os.makedirs(jdir, mode=0o750, exist_ok=True)
    out({"ok": True, "job_id": args[0], "path": jdir})


def cmd_put(cfg, args):
    if len(args) != 2:
        raise HpcError("usage: put <job_id> <file>")
    jdir = job_dir(cfg, args[0])
    name = args[1]
    if not UPLOAD_NAMES.match(name):
        raise HpcError("file name not allowed: " + name)
    if os.path.exists(os.path.join(jdir, "state.json")):
        raise HpcError("job already submitted")
    limit = max(1, ival(cfg, "MAX_UPLOAD_MB")) * 1024 * 1024
    dest = os.path.join(jdir, name)
    tmp = dest + ".part"
    size = 0
    with open(tmp, "wb") as fh:
        src = sys.stdin.buffer if hasattr(sys.stdin, "buffer") else sys.stdin
        while True:
            chunk = src.read(1 << 16)
            if not chunk:
                break
            size += len(chunk)
            if size > limit:
                fh.close()
                os.remove(tmp)
                raise HpcError("file too large (limit %d MB)" % (limit // 1048576))
            fh.write(chunk)
    os.replace(tmp, dest)
    out({"ok": True, "file": name, "bytes": size})


def validate_spec(cfg, jdir):
    spec = read_json(os.path.join(jdir, "job.json"))
    if not isinstance(spec, dict):
        raise HpcError("job.json missing or invalid")
    mode = spec.get("mode")
    if mode == "forward":
        f = spec.get("forward") or {}
        if not os.path.isfile(os.path.join(jdir, f.get("receptor", "receptor.pdb"))):
            raise HpcError("receptor.pdb missing")
        src = f.get("source", "database")
        if src in ("single", "custom"):
            mf = f.get("molecule_file", "")
            if not UPLOAD_NAMES.match(mf or "") or not os.path.isfile(os.path.join(jdir, mf)):
                raise HpcError("molecule file missing")
        elif src == "database":
            info = None
            for e in catalog(cfg, "databases"):
                if e["key"] == str(f.get("database", "")).lower():
                    info = dataset_info(cfg, e, "databases")
            if not info:
                raise HpcError("unknown compound library")
            if not info["available"]:
                raise HpcError("library not installed: " + info["label"])
        else:
            raise HpcError("invalid source")
    elif mode == "reverse":
        r = spec.get("reverse") or {}
        qf = r.get("query_file", "")
        if not UPLOAD_NAMES.match(qf or "") or not os.path.isfile(os.path.join(jdir, qf)):
            raise HpcError("query molecule missing")
        info = None
        for e in catalog(cfg, "libraries"):
            if e["key"] == str(r.get("library", "homo")).lower():
                info = dataset_info(cfg, e, "libraries")
        if not info:
            raise HpcError("unknown target library")
        if not info["available"]:
            raise HpcError("target library not installed: " + info["label"])
    else:
        raise HpcError("invalid mode")
    return spec


def cmd_submit(cfg, args):
    if len(args) != 1:
        raise HpcError("usage: submit <job_id>")
    jdir = job_dir(cfg, args[0])
    if os.path.exists(os.path.join(jdir, "state.json")):
        st = read_json(os.path.join(jdir, "state.json"), {}) or {}
        out({"ok": True, "state": st.get("state"), "already_submitted": True,
             "scheduler_id": st.get("scheduler_id")})
    spec = validate_spec(cfg, jdir)
    profile = job_profile(cfg, spec)
    opts = spec.setdefault("options", {})
    if profile == "large" and ival(cfg, "RESULT_LIMIT_LARGE") and not opts.get("top"):
        opts["top"] = ival(cfg, "RESULT_LIMIT_LARGE")
    opts.setdefault("preview_rows", ival(cfg, "PREVIEW_ROWS"))
    if spec.get("mode") == "forward" and (spec.get("forward") or {}).get("source") == "custom":
        opts.setdefault("max_molecules", ival(cfg, "MAX_CUSTOM_MOLECULES"))
    write_json(os.path.join(jdir, "job.json"), spec)
    state = {"state": "QUEUED", "queued_at": now_iso(), "updated_at": now_iso(),
             "scheduler": cfg["SCHEDULER"], "profile": profile}
    write_json(os.path.join(jdir, "state.json"), state)
    open(os.path.join(jdir, "QUEUED"), "a").close()
    try:
        if cfg["SCHEDULER"] == "slurm":
            sid = submit_slurm(cfg, jdir, profile)
        elif cfg["SCHEDULER"] == "pbs":
            sid = submit_pbs(cfg, jdir, profile)
        else:
            sid = submit_local(cfg, jdir, profile)
    except HpcError as exc:
        mark_failed(jdir, "The job could not be submitted to the cluster: %s" % exc)
        raise
    st = read_json(os.path.join(jdir, "state.json"), {}) or {}
    st["scheduler_id"] = sid
    write_json(os.path.join(jdir, "state.json"), st)
    out({"ok": True, "job_id": args[0], "state": st.get("state", "QUEUED"), "scheduler": cfg["SCHEDULER"],
         "scheduler_id": sid, "profile": profile})


def legacy_progress(jdir):
    """Aggregate the 32-partition progress files written by the original HPC script."""
    total = cur = active = 0
    for i in range(1, 65):
        pj = read_json(os.path.join(jdir, str(i), "progress.json"))
        if pj:
            active += 1
            cur += int(pj.get("current", 0) or 0)
            total += int(pj.get("total", 0) or 0)
    if not active:
        return None
    pct = min(100, int(100 * cur / total)) if total else 0
    return {"stage": "scoring", "stage_label": "Parallel screening", "percent": pct, "current": cur,
            "total": total, "message": "Running on %d partitions" % active}


def status_payload(cfg, job_id, check_scheduler=True):
    jdir = job_dir(cfg, job_id)
    st = read_json(os.path.join(jdir, "state.json"))
    legacy = st is None
    if legacy:
        st = {}
        for flag in ("COMPLETED", "FAILED", "CANCELLED", "RUNNING", "QUEUED"):
            if os.path.exists(os.path.join(jdir, flag)):
                st["state"] = flag
                break
        st.setdefault("state", "UNKNOWN")
    state = st.get("state", "UNKNOWN")
    prog = read_json(os.path.join(jdir, "progress.json")) or (legacy_progress(jdir) if legacy else None) or {}

    sched_state, reason = None, ""
    if check_scheduler and state in ("QUEUED", "RUNNING") and st.get("scheduler_id"):
        sched_state, reason = scheduler_state(cfg, st.get("scheduler", cfg["SCHEDULER"]), st.get("scheduler_id"))
        gone = sched_state is None or str(sched_state).startswith("ENDED:")
        if gone:
            time.sleep(0.5)  # the pipeline may be writing its final state right now
            st2 = read_json(os.path.join(jdir, "state.json"), {}) or {}
            if st2.get("state") in TERMINAL:
                st, state = st2, st2["state"]
            else:
                ended = str(sched_state or "").replace("ENDED:", "")
                quiet_for = time.time() - os.path.getmtime(os.path.join(jdir, "state.json"))
                # squeue/qstat can briefly miss a job right after submission; require a quiet period
                if ended or cfg["SCHEDULER"] == "local" or quiet_for > 60:
                    if os.path.exists(os.path.join(jdir, "CANCEL_REQUESTED")):
                        mark_cancelled(jdir)
                    else:
                        mark_failed(jdir, SLURM_END_MESSAGES.get(
                            ended, "The job stopped unexpectedly on the cluster."))
                    st = read_json(os.path.join(jdir, "state.json"), {}) or {}
                    state = st.get("state", "FAILED")

    payload = {
        "ok": True,
        "job_id": job_id,
        "state": state,
        "legacy": legacy,
        "stage": prog.get("stage"),
        "stage_label": prog.get("stage_label"),
        "percent": prog.get("percent"),
        "message": prog.get("message"),
        "current": prog.get("current"),
        "total": prog.get("total"),
        "scored": prog.get("scored"),
        "eta_seconds": prog.get("eta_seconds"),
        "rate_per_s": prog.get("rate_per_s"),
        "queued_at": st.get("queued_at"),
        "started_at": st.get("started_at"),
        "finished_at": st.get("finished_at"),
        "elapsed_seconds": st.get("elapsed_seconds"),
        "scheduler": st.get("scheduler"),
        "scheduler_id": st.get("scheduler_id"),
        "scheduler_state": sched_state,
        "scheduler_reason": reason or None,
        "error": None,
        "results": None,
    }
    if state == "QUEUED" and not prog:
        payload["stage_label"] = "Waiting in the cluster queue"
        if reason and reason not in ("None", "(null)"):
            payload["message"] = "Scheduler reason: " + reason
    if state in ("FAILED", "CANCELLED"):
        text = tail(os.path.join(jdir, "error.log"), 200, 32768)
        msg, suggestion, details = (st.get("message") or "The job failed."), None, ""
        if text:
            head, _, rest = text.partition("\n--- Details ---\n")
            lines = [ln for ln in head.splitlines() if ln.strip()]
            if lines:
                msg = lines[0].strip()
            for ln in lines:
                if ln.startswith("Suggestion:"):
                    suggestion = ln.split(":", 1)[1].strip()
            details = rest[-4000:]
        payload["error"] = {"message": msg, "suggestion": suggestion, "details": details,
                            "exit_code": st.get("exit_code"), "stage": st.get("stage")}
    res = {}
    for name, key in (("results.csv", "full"), ("results.csv.gz", "gz"), ("results_preview.csv", "preview"),
                      ("summary.json", "summary")):
        path = os.path.join(jdir, name)
        if os.path.isfile(path):
            res[key] = os.path.getsize(path)
    if res:
        summ = read_json(os.path.join(jdir, "summary.json"), {}) or {}
        res["rows"] = summ.get("rows_written")
        payload["results"] = res
    return payload


def cmd_status(cfg, args):
    if len(args) < 1:
        raise HpcError("usage: status <job_id>")
    out(status_payload(cfg, args[0], check_scheduler="--no-scheduler" not in args))


def cmd_get(cfg, args):
    if len(args) != 2:
        raise HpcError("usage: get <job_id> <file>")
    jdir = job_dir(cfg, args[0])
    name = args[1]
    if name not in DOWNLOAD_NAMES:
        raise HpcError("file not available: " + name)
    path = os.path.join(jdir, name)
    if not os.path.isfile(path):
        raise HpcError("file not found: " + name)
    dst = sys.stdout.buffer if hasattr(sys.stdout, "buffer") else sys.stdout
    with open(path, "rb") as fh:
        shutil.copyfileobj(fh, dst, 1 << 20)
    dst.flush()
    sys.exit(0)


def cmd_cancel(cfg, args):
    if len(args) != 1:
        raise HpcError("usage: cancel <job_id>")
    jdir = job_dir(cfg, args[0])
    st = read_json(os.path.join(jdir, "state.json"), {}) or {}
    state = st.get("state")
    if state in TERMINAL:
        out({"ok": True, "state": state, "note": "job already finished"})
    open(os.path.join(jdir, "CANCEL_REQUESTED"), "a").close()
    sched, sid = st.get("scheduler", cfg["SCHEDULER"]), st.get("scheduler_id", "")
    if sched == "slurm" and sid:
        run(["scancel", "--signal=TERM", "--full", sid])
        if state == "QUEUED":
            run(["scancel", sid])
            mark_cancelled(jdir)
    elif sched == "pbs" and sid:
        run(["qdel", sid])
        if state == "QUEUED":
            mark_cancelled(jdir)
    elif sched == "local" and sid.startswith("local-"):
        pid = int(sid.split("-", 1)[1])
        try:
            os.killpg(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            pass
        if state == "QUEUED":
            mark_cancelled(jdir)
    else:
        mark_cancelled(jdir)
    for _ in range(10):
        time.sleep(0.5)
        st = read_json(os.path.join(jdir, "state.json"), {}) or {}
        if st.get("state") in TERMINAL:
            break
    out({"ok": True, "state": st.get("state", "CANCELLING")})


def cmd_list(cfg, args):
    limit = 30
    if "--limit" in args:
        limit = int(args[args.index("--limit") + 1])
    rows = []
    base = cfg["JOBS_DIR"]
    if os.path.isdir(base):
        names = [n for n in os.listdir(base) if JOB_ID_RE.match(n) and os.path.isdir(os.path.join(base, n))]
        names.sort(key=lambda n: os.path.getmtime(os.path.join(base, n)), reverse=True)
        for n in names[:limit]:
            st = read_json(os.path.join(base, n, "state.json"), {}) or {}
            spec = read_json(os.path.join(base, n, "job.json"), {}) or {}
            rows.append({"job_id": n, "state": st.get("state", "?"), "mode": spec.get("mode"),
                         "queued_at": st.get("queued_at"), "finished_at": st.get("finished_at")})
    out({"ok": True, "jobs": rows})


def cmd_cleanup(cfg, args):
    days = ival(cfg, "JOB_RETENTION_DAYS")
    if "--days" in args:
        days = int(args[args.index("--days") + 1])
    dry = "--dry-run" in args
    if days <= 0:
        out({"ok": True, "removed": [], "note": "retention disabled"})
    cutoff = time.time() - days * 86400
    removed = []
    base = cfg["JOBS_DIR"]
    for n in os.listdir(base) if os.path.isdir(base) else []:
        path = os.path.join(base, n)
        if not JOB_ID_RE.match(n) or not os.path.isdir(path):
            continue
        st = read_json(os.path.join(path, "state.json"), {}) or {}
        terminal = st.get("state") in TERMINAL or any(
            os.path.exists(os.path.join(path, f)) for f in ("COMPLETED", "FAILED", "CANCELLED"))
        mtime = os.path.getmtime(path)
        if (terminal and mtime < cutoff) or mtime < cutoff - days * 86400:
            removed.append(n)
            if not dry:
                shutil.rmtree(path, ignore_errors=True)
    out({"ok": True, "removed": removed, "dry_run": dry, "days": days})


def cmd_selftest(cfg, args):
    simulate = "--real" not in args
    job_id = "selftest-" + datetime.now().strftime("%Y%m%d%H%M%S")
    jdir = job_dir(cfg, job_id, must_exist=False)
    os.makedirs(jdir)
    example = os.path.join(cfg["SEARCH_ML_HOME"], "example", "9kte.pdb")
    shutil.copyfile(example, os.path.join(jdir, "receptor.pdb"))
    spec = {"version": 2, "mode": "forward",
            "forward": {"receptor": "receptor.pdb", "ligand": {"resname": "9UM", "chain": "A", "resseq": "413"},
                        "source": "database", "database": "fda", "filters": {}},
            "options": {"simulate": simulate, "simulate_seconds": 4}}
    write_json(os.path.join(jdir, "job.json"), spec)
    t0 = time.time()
    rc, o, e = run([cfg["PYTHON"], os.path.join(cfg["SEARCH_ML_HOME"], "scripts", "pipeline.py"), "run", jdir,
                    "--quiet"], timeout=3600)
    status = status_payload(cfg, job_id, check_scheduler=False)
    top = []
    prev = os.path.join(jdir, "results_preview.csv")
    if os.path.isfile(prev):
        with open(prev) as fh:
            top = [ln.strip() for ln in fh.readlines()[:6]]
    out({"ok": rc == 0, "job_id": job_id, "simulated": simulate, "exit_code": rc,
         "seconds": round(time.time() - t0, 1), "state": status["state"], "error": status["error"],
         "top_rows": top, "stderr": e[-800:]}, 0 if rc == 0 else 1)


# --------------------------------------------------------------------------
# SSH forced-command gateway
# --------------------------------------------------------------------------

GATEWAY_COMMANDS = {"ping", "prepare", "put", "submit", "status", "get", "cancel", "cleanup", "version"}


def gateway(cfg, original):
    """Validate an SSH_ORIGINAL_COMMAND and run it (only searchml-hpc commands and scp into JOBS_DIR)."""
    try:
        tokens = shlex.split(original or "")
    except ValueError:
        fail("rejected: unparsable command")
    if not tokens:
        fail("rejected: interactive logins are disabled for this key")
    if tokens[0] == "scp":
        flags = [t for t in tokens[1:] if t.startswith("-")]
        paths = [t for t in tokens[1:] if not t.startswith("-")]
        allowed = {"-t", "-f", "-p", "-d", "-q", "--"}
        if not set(flags) <= allowed or len(paths) != 1 or ("-t" in flags) == ("-f" in flags):
            fail("rejected: scp options not allowed")
        target = os.path.realpath(paths[0])
        jobs = os.path.realpath(cfg["JOBS_DIR"])
        if not target.startswith(jobs + os.sep):
            fail("rejected: path outside the jobs directory")
        rel = os.path.relpath(target, jobs).split(os.sep)
        if not JOB_ID_RE.match(rel[0]):
            fail("rejected: bad job directory")
        if "-t" in flags and len(rel) == 2 and not UPLOAD_NAMES.match(rel[1]):
            fail("rejected: upload file name not allowed")
        if "-f" in flags and (len(rel) != 2 or rel[1] not in DOWNLOAD_NAMES):
            fail("rejected: download not allowed")
        os.execvp("scp", ["scp"] + flags + [target])
    if os.path.basename(tokens[0]) in ("searchml-hpc", "searchml_hpc.py"):
        tokens = tokens[1:]
    if not tokens or tokens[0] not in GATEWAY_COMMANDS:
        fail("rejected: command not allowed")
    dispatch(cfg, tokens)


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------

def dispatch(cfg, argv):
    cmd, args = argv[0], argv[1:]
    table = {
        "ping": cmd_ping, "prepare": cmd_prepare, "put": cmd_put, "submit": cmd_submit,
        "status": cmd_status, "get": cmd_get, "cancel": cmd_cancel, "list": cmd_list,
        "cleanup": cmd_cleanup, "selftest": cmd_selftest,
    }
    if cmd == "version":
        out({"ok": True, "version": VERSION})
    if cmd == "_local_worker":
        sys.exit(local_worker(cfg, args[0], args[1] if len(args) > 1 else "small"))
    if cmd not in table:
        fail("unknown command '%s' (try: %s)" % (cmd, ", ".join(sorted(table))))
    try:
        table[cmd](cfg, args)
    except HpcError as exc:
        fail(str(exc))
    except SystemExit:
        raise
    except Exception as exc:  # pragma: no cover
        fail("internal error: %s: %s" % (type(exc).__name__, exc))


def main():
    argv = sys.argv[1:]
    cfg = load_config()
    if argv and argv[0] == "--gateway":
        gateway(cfg, argv[1] if len(argv) > 1 else os.environ.get("SSH_ORIGINAL_COMMAND", ""))
        return
    if not argv or argv[0] in ("-h", "--help", "help"):
        sys.stdout.write(__doc__)
        return
    dispatch(cfg, argv)


if __name__ == "__main__":
    main()
