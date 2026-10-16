#!/usr/bin/env python
"""
SEARCH-ML pipeline orchestrator.

One code path for every way of running SEARCH-ML:

* ``master_screening.sh`` / ``master_screening_reverse.sh`` (command line)
* the HPC job runner used by the web portal (``hpc/bin/searchml-hpc``)

Sub-commands
------------
    pipeline.py forward --pdb 9kte.pdb --ligand 9UM --mode DB
    pipeline.py forward --pdb 9kte.pdb --ligand 9UM:A:413 --mode SINGLE --molecules aspirin.smi
    pipeline.py forward --pdb 9kte.pdb --ligand 9UM --mode CUSTOM --molecules library.sdf --no-filters
    pipeline.py reverse --query query_ligand.sdf --library HOMO
    pipeline.py reverse --smiles "CC(=O)Oc1ccccc1C(=O)O"
    pipeline.py run <job_dir>                 # execute a prepared job directory (job.json)
    pipeline.py check                         # verify the environment (RDKit, Biopython, tleap, models)
    pipeline.py catalog                       # list compound libraries / target libraries

Job directory layout
--------------------
    job.json                 job specification (inputs, mode, library, filters)
    state.json               QUEUED | RUNNING | COMPLETED | FAILED | CANCELLED (+ timestamps)
    RUNNING / COMPLETED / FAILED   flag files (kept for backwards compatibility)
    progress.json            live progress (stage, percent, rows scored, ETA)
    job.log                  full log of every step
    error.log                user-facing error message + technical details (on failure)
    results.csv              ranked predictions (strongest binders first)
    results_preview.csv      top rows of results.csv (for the web portal)
    summary.json             counts, filters, affinity histogram over all scored rows
    PROT.pdb, PROT_LIG.pdb, leap.log, protein.csv, dataset.csv ...   intermediates

Author: Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi
"""

from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import json
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import time
from collections import deque
from typing import Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)

__version__ = "2.0.0"

# --------------------------------------------------------------------------
# Exit codes (1-6 identical to the original master scripts)
# --------------------------------------------------------------------------
EXIT_HOME = 1
EXIT_INPUT_MISSING = 2
EXIT_PROTEIN_SHORT = 3
EXIT_SINGLE_MISSING = 4
EXIT_CUSTOM_MISSING = 5
EXIT_BAD_MODE = 6
EXIT_MODELS = 7
EXIT_FEATURES = 8
EXIT_MEMORY = 9
EXIT_PROTEIN_FEATURES = 10
EXIT_POCKET_FEATURES = 11
EXIT_MERGE = 12
EXIT_LIGAND_NOT_FOUND = 13
EXIT_TLEAP = 14
EXIT_MOLECULE = 15
EXIT_NO_VALID_MOLECULES = 16
EXIT_ENVIRONMENT = 17
EXIT_LIBRARY_MISSING = 18
EXIT_INTERNAL = 20
EXIT_CANCELLED = 130

FILTERABLE = [
    "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds",
    "NumAromaticRings", "NumHeteroatoms", "HeavyAtomCount", "NumAtoms", "MolMR", "NumCarbons",
    "FractionCSP3", "BertzCT", "LigVol", "RadiusOfGyration",
]

SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,80}$")


class PipelineError(Exception):
    def __init__(self, message: str, exit_code: int = EXIT_INTERNAL, details: str = "", hint: str = ""):
        super().__init__(message)
        self.exit_code = exit_code
        self.details = details
        self.hint = hint


class Cancelled(Exception):
    pass


# --------------------------------------------------------------------------
# Paths & catalogs
# --------------------------------------------------------------------------

def search_ml_home() -> str:
    home = os.environ.get("SEARCH_ML_HOME")
    if home and os.path.isdir(home):
        return os.path.realpath(home)
    return os.path.dirname(HERE)  # scripts/.. (repository root)


def datasets_dir() -> str:
    return os.environ.get("SEARCHML_DATASETS") or os.path.join(search_ml_home(), "datasets")


def models_dir() -> str:
    return os.environ.get("SEARCHML_MODELS") or os.path.join(search_ml_home(), "models")


def parameters_dir() -> str:
    return os.path.join(search_ml_home(), "parameters")


def _catalog_path(kind: str) -> str:
    env = os.environ.get("SEARCHML_DATABASES_TSV" if kind == "databases" else "SEARCHML_LIBRARIES_TSV")
    if env:
        return env
    return os.path.join(search_ml_home(), "config", f"{kind}.tsv")


def load_catalog(kind: str) -> Dict[str, dict]:
    """Read config/databases.tsv or config/libraries.tsv -> {key: entry}."""
    path = _catalog_path(kind)
    entries: Dict[str, dict] = {}
    if not os.path.isfile(path):
        return entries
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            key = cols[0].strip().lower()
            entries[key] = {
                "key": key,
                "aliases": [a.strip() for a in cols[1].split(",") if a.strip()],
                "label": cols[2].strip(),
                "path": cols[3].strip(),
                "size": (cols[4].strip() if len(cols) > 4 else "small") or "small",
            }
    return entries


def resolve_catalog_key(kind: str, name: str) -> Optional[dict]:
    cat = load_catalog(kind)
    if not name:
        return None
    low = name.strip().lower()
    if low in cat:
        return cat[low]
    for e in cat.values():
        if any(low == a.lower() for a in e["aliases"]):
            return e
    return None


def resolve_dataset_inputs(entry: dict, kind: str) -> Tuple[List[str], Optional[int]]:
    """Return ([files or dir], known_row_count) for a catalog entry."""
    path = entry["path"]
    if not os.path.isabs(path):
        path = os.path.join(datasets_dir(), path)
    main_name = "dataset.csv" if kind == "databases" else "protein.csv"
    count = None
    if os.path.isdir(path):
        count_file = os.path.join(path, "count.txt")
        if os.path.isfile(count_file):
            try:
                count = int(open(count_file).read().split()[0])
            except (ValueError, IndexError):
                count = None
        for candidate in (main_name, main_name + ".gz"):
            if os.path.isfile(os.path.join(path, candidate)):
                return [os.path.join(path, candidate)], count
        parts = os.path.join(path, "parts")
        if os.path.isdir(parts) and any(f.endswith((".csv", ".csv.gz")) for f in os.listdir(parts)):
            return [parts], count
        if any(f.endswith((".csv", ".csv.gz")) for f in os.listdir(path)):
            return [path], count
        raise PipelineError(f"Library '{entry['label']}' has no data files in {path}.", EXIT_LIBRARY_MISSING)
    if os.path.isfile(path):
        return [path], count
    raise PipelineError(
        f"Library '{entry['label']}' is not installed on this server (expected at {path}).",
        EXIT_LIBRARY_MISSING,
        hint="Ask the administrator to install the dataset or edit config/%s.tsv." % kind)


# --------------------------------------------------------------------------
# Small utilities
# --------------------------------------------------------------------------

def now_iso() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


def atomic_write_json(path: str, payload: dict) -> None:
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def read_json(path: str, default=None):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return default


def touch(path: str) -> None:
    with open(path, "a"):
        os.utime(path, None)


def rm_f(path: str) -> None:
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


def tail_text(path: str, n: int = 25) -> str:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return "".join(deque(fh, maxlen=n))
    except OSError:
        return ""


def detect_cpus() -> int:
    for var in ("SEARCHML_CPUS", "SLURM_CPUS_PER_TASK", "NCPUS", "PBS_NUM_PPN", "LSB_DJOB_NUMPROC"):
        v = os.environ.get(var)
        if v and v.isdigit() and int(v) > 0:
            return int(v)
    try:
        return max(1, len(os.sched_getaffinity(0)))  # type: ignore[attr-defined]
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def find_tleap() -> Optional[str]:
    path = shutil.which("tleap")
    if path:
        return path
    for base in filter(None, [os.environ.get("AMBERHOME"), sys.prefix]):
        cand = os.path.join(base, "bin", "tleap")
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


def check_environment(mode: str, verbose: bool = False) -> List[str]:
    """Return a list of problems (empty = OK)."""
    problems = []
    try:
        import rdkit  # noqa: F401
    except ImportError:
        problems.append("RDKit is not importable (install the search-ml conda environment).")
    for mod, label in (("pandas", "pandas"), ("numpy", "numpy"), ("joblib", "joblib"),
                       ("sklearn", "scikit-learn")):
        try:
            __import__(mod)
        except ImportError:
            problems.append(f"{label} is not importable.")
    if mode == "forward":
        try:
            import Bio  # noqa: F401
        except ImportError:
            problems.append("Biopython is not importable.")
        try:
            import scipy  # noqa: F401
        except ImportError:
            problems.append("SciPy is not importable.")
        if not find_tleap():
            problems.append("AmberTools 'tleap' was not found on PATH (activate the search-ml env).")
    for mod in ("xgboost", "catboost", "lightgbm"):
        try:
            __import__(mod)
        except ImportError:
            problems.append(f"{mod} is not importable.")
    mdir = models_dir()
    import scoring
    needed = [scoring.SCALER_FILE, scoring.META_FILE, *scoring.MODEL_FILES.values()]
    missing = [f for f in needed if not os.path.isfile(os.path.join(mdir, f))]
    if missing:
        problems.append(f"Missing model files in {mdir}: {', '.join(missing)}")
    if verbose:
        print(f"SEARCH_ML_HOME : {search_ml_home()}")
        print(f"Python         : {sys.executable} ({sys.version.split()[0]})")
        print(f"tleap          : {find_tleap() or 'NOT FOUND'}")
        print(f"Models         : {mdir}")
        print(f"Datasets       : {datasets_dir()}")
    return problems


# --------------------------------------------------------------------------
# Job runner
# --------------------------------------------------------------------------

class JobRunner:
    """Executes one job directory described by job.json."""

    def __init__(self, job_dir: str, workers: Optional[int] = None, threads: Optional[int] = None,
                 echo: bool = True, simulate: bool = False):
        self.dir = os.path.realpath(job_dir)
        self.spec_path = os.path.join(self.dir, "job.json")
        self.echo = echo
        self.log_fh = None
        self.child: Optional[subprocess.Popen] = None
        self.cancel_requested = False
        self.stage = "starting"
        self.spec = read_json(self.spec_path)
        if not isinstance(self.spec, dict):
            raise PipelineError(f"{self.spec_path} is missing or not valid JSON.", EXIT_INPUT_MISSING)
        opts = self.spec.get("options") or {}
        self.simulate = bool(simulate or opts.get("simulate") or os.environ.get("SEARCHML_SIMULATE") == "1")
        cpus = detect_cpus()
        self.cpus = cpus
        self.workers = workers or int(opts.get("workers") or 0) or None
        self.threads = threads

    # -- paths / logging ---------------------------------------------------
    def p(self, name: str) -> str:
        return os.path.join(self.dir, name)

    def log(self, msg: str) -> None:
        line = f"[{_dt.datetime.now().strftime('%H:%M:%S')}] {msg}"
        if self.log_fh:
            self.log_fh.write(line + "\n")
            self.log_fh.flush()
        if self.echo:
            print(line, flush=True)

    # -- state -------------------------------------------------------------
    def set_state(self, state: str, **extra) -> None:
        current = read_json(self.p("state.json"), {}) or {}
        current.update(extra)
        current["state"] = state
        current["updated_at"] = now_iso()
        atomic_write_json(self.p("state.json"), current)
        for flag in ("QUEUED", "RUNNING", "COMPLETED", "FAILED", "CANCELLED"):
            if flag != state:
                rm_f(self.p(flag))
        if state in ("RUNNING", "COMPLETED", "FAILED", "CANCELLED"):
            touch(self.p(state))

    def progress(self, stage: str, label: str, percent: float, message: str = "") -> None:
        self.stage = stage
        atomic_write_json(self.p("progress.json"), {
            "stage": stage,
            "stage_label": label,
            "percent": round(percent, 1),
            "message": message or label,
            "updated_at": int(time.time()),
        })
        self.log(f"== {label} ({percent:.0f}%)")

    # -- subprocess --------------------------------------------------------
    def run_step(self, args: List[str], error_message: str, exit_code: int,
                 timeout: Optional[int] = None, hint: str = "", cwd: Optional[str] = None,
                 env: Optional[dict] = None, code_map: Optional[Dict[int, int]] = None) -> str:
        """Run a command, streaming its output to job.log. Raises PipelineError on failure."""
        if self.cancel_requested:
            raise Cancelled()
        self.log("$ " + " ".join(_shell_quote(a) for a in args))
        tail: deque = deque(maxlen=40)
        try:
            self.child = subprocess.Popen(
                args, cwd=cwd or self.dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                env=env or os.environ.copy(), text=True, errors="replace", bufsize=1)
        except FileNotFoundError as exc:
            raise PipelineError(error_message, exit_code, details=str(exc), hint=hint)
        start = time.time()
        assert self.child.stdout is not None
        for line in self.child.stdout:
            tail.append(line)
            if self.log_fh:
                self.log_fh.write(line)
            if self.echo:
                sys.stdout.write(line)
            if timeout and time.time() - start > timeout:
                self.child.kill()
                raise PipelineError(f"{error_message} (step exceeded {timeout // 60} min)", exit_code,
                                    details="".join(tail), hint=hint)
        rc = self.child.wait()
        self.child = None
        if self.log_fh:
            self.log_fh.flush()
        if self.cancel_requested:
            raise Cancelled()
        if rc != 0:
            code = (code_map or {}).get(rc, exit_code)
            reason = _last_error_line(tail)
            if rc in (-9, 137) or (rc == EXIT_MEMORY and code == exit_code):
                code, reason = EXIT_MEMORY, "The process ran out of memory or was killed."
            message = f"{error_message} {reason}".strip() if reason else error_message
            raise PipelineError(message, code, details="".join(tail), hint=hint)
        return "".join(tail)

    def py(self, script: str) -> str:
        return os.path.join(HERE, script)

    # -- spec helpers ------------------------------------------------------
    def _input(self, name: Optional[str], code: int, what: str) -> str:
        if not name or not SAFE_NAME.match(name) or ".." in name:
            raise PipelineError(f"Invalid or missing {what} file name in job.json.", code)
        path = self.p(name)
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise PipelineError(f"The {what} file '{name}' is missing or empty.", code)
        return path

    def _write_filters(self, filters) -> Optional[str]:
        """filters: None -> default parameter file, {} -> no filters, dict -> those ranges."""
        dest = self.p("parameter.txt")
        if filters is None:
            default = os.path.join(parameters_dir(), "parameter.txt")
            if os.path.isfile(default):
                shutil.copyfile(default, dest)
                return dest
            return None
        if not isinstance(filters, dict):
            raise PipelineError("'filters' in job.json must be an object.", EXIT_BAD_MODE)
        lines = []
        for name, rng in filters.items():
            if name not in FILTERABLE:
                raise PipelineError(f"Unknown filter property '{name}'.", EXIT_BAD_MODE)
            try:
                lo, hi = float(rng[0]), float(rng[1])
            except (TypeError, ValueError, IndexError):
                raise PipelineError(f"Filter '{name}' needs [min, max] numbers.", EXIT_BAD_MODE)
            if lo > hi:
                lo, hi = hi, lo
            lines.append(f"{name}: {lo:g},{hi:g}")
        with open(dest, "w") as fh:
            fh.write("\n".join(lines) + ("\n" if lines else ""))
        return dest if lines else None

    # -- main --------------------------------------------------------------
    def execute(self) -> int:
        os.makedirs(self.dir, exist_ok=True)
        self.log_fh = open(self.p("job.log"), "a", encoding="utf-8")
        rm_f(self.p("error.log"))
        rm_f(self.p("CANCEL_REQUESTED"))
        signal.signal(signal.SIGTERM, self._on_signal)
        signal.signal(signal.SIGINT, self._on_signal)
        t0 = time.time()
        mode = (self.spec.get("mode") or "").lower()
        self.set_state("RUNNING", started_at=now_iso(), pid=os.getpid(),
                       host=os.uname().nodename, mode=mode)
        self.log(f"SEARCH-ML pipeline {__version__} | job {os.path.basename(self.dir)} | mode {mode}"
                 f"{' | SIMULATION' if self.simulate else ''}")
        self.log(f"Host {os.uname().nodename} | CPUs {self.cpus} | Python {sys.version.split()[0]}")
        try:
            if mode not in ("forward", "reverse"):
                raise PipelineError(f"Unknown mode '{mode}' (expected forward or reverse).", EXIT_BAD_MODE)
            if not self.simulate:
                problems = check_environment(mode)
                if problems:
                    raise PipelineError("The compute environment is not set up correctly.", EXIT_ENVIRONMENT,
                                        details="\n".join(problems),
                                        hint="Run 'python scripts/pipeline.py check' on the server.")
            if self.simulate:
                summary = self._simulate(mode)
            elif mode == "forward":
                summary = self._forward()
            else:
                summary = self._reverse()
            self._finalize(summary)
            elapsed = round(time.time() - t0, 1)
            self.progress("done", "Completed", 100.0, "Results are ready")
            self.set_state("COMPLETED", finished_at=now_iso(), elapsed_seconds=elapsed, exit_code=0,
                           rows=summary.get("rows_written"))
            self.log(f"Job completed successfully in {elapsed:.1f} s.")
            return 0
        except Cancelled:
            self.log("Job cancelled.")
            self._write_error("The job was cancelled.", "", "", EXIT_CANCELLED)
            self.set_state("CANCELLED", finished_at=now_iso(), exit_code=EXIT_CANCELLED, stage=self.stage)
            return EXIT_CANCELLED
        except PipelineError as exc:
            self.log(f"ERROR: {exc}")
            self._write_error(str(exc), exc.details, exc.hint, exc.exit_code)
            self.set_state("FAILED", finished_at=now_iso(), exit_code=exc.exit_code, stage=self.stage,
                           message=str(exc))
            return exc.exit_code
        except MemoryError:
            msg = "The job ran out of memory."
            self._write_error(msg, "", "Try a smaller library or ask the administrator for more memory.",
                              EXIT_MEMORY)
            self.set_state("FAILED", finished_at=now_iso(), exit_code=EXIT_MEMORY, stage=self.stage, message=msg)
            return EXIT_MEMORY
        except Exception as exc:  # unexpected
            import traceback
            tb = traceback.format_exc()
            self.log(tb)
            msg = "Internal pipeline error."
            self._write_error(msg, tb, "Please report this job ID to the SEARCH-ML team.", EXIT_INTERNAL)
            self.set_state("FAILED", finished_at=now_iso(), exit_code=EXIT_INTERNAL, stage=self.stage,
                           message=f"{msg} {exc}")
            return EXIT_INTERNAL
        finally:
            if self.log_fh:
                self.log_fh.close()
                self.log_fh = None

    def _on_signal(self, signum, _frame):
        self.cancel_requested = True
        if os.path.exists(self.p("CANCEL_REQUESTED")) or signum == signal.SIGINT:
            reason = "cancel"
        else:
            reason = "terminated"
        if self.child and self.child.poll() is None:
            try:
                self.child.terminate()
            except Exception:
                pass
        if reason == "terminated":
            # scheduler time limit / node shutdown: report as a failure, not a user cancel
            raise PipelineError("The job was stopped by the cluster scheduler (time limit or node shutdown).",
                                EXIT_CANCELLED, hint="Large libraries may need a longer wall-time limit.")
        raise Cancelled()

    def _write_error(self, message: str, details: str, hint: str, code: int) -> None:
        with open(self.p("error.log"), "w", encoding="utf-8") as fh:
            fh.write(message.strip() + "\n")
            if hint:
                fh.write("\nSuggestion: " + hint.strip() + "\n")
            fh.write(f"\nStage: {self.stage}\nExit code: {code}\n")
            if details:
                fh.write("\n--- Details ---\n" + details.strip()[-6000:] + "\n")

    # -- forward screening -------------------------------------------------
    def _forward(self) -> dict:
        f = self.spec.get("forward") or {}
        receptor = self._input(f.get("receptor", "receptor.pdb"), EXIT_INPUT_MISSING, "protein structure")
        lig = f.get("ligand") or {}
        if isinstance(lig, dict):
            selector = f"{lig.get('resname', '')}:{lig.get('chain') or '_'}:{lig.get('resseq', '')}{lig.get('icode', '')}"
            if not lig.get("resseq"):
                selector = str(lig.get("resname", ""))
        else:
            selector = str(lig)
        source = (f.get("source") or "database").lower()
        py = sys.executable

        # resolve the library early so a missing dataset fails fast
        mol_inputs: List[str] = []
        known_rows = None
        db_entry = None
        if source == "database":
            db_entry = resolve_catalog_key("databases", f.get("database", ""))
            if not db_entry:
                raise PipelineError(f"Unknown compound library '{f.get('database')}'.", EXIT_BAD_MODE)
            mol_inputs, known_rows = resolve_dataset_inputs(db_entry, "databases")
        elif source == "single":
            mol_file = self._input(f.get("molecule_file"), EXIT_SINGLE_MISSING, "molecule")
        elif source == "custom":
            mol_file = self._input(f.get("molecule_file"), EXIT_CUSTOM_MISSING, "molecule library")
        else:
            raise PipelineError(f"Invalid screening source '{source}'.", EXIT_BAD_MODE)

        params = self._write_filters(f.get("filters") if "filters" in f else
                                     (None if source == "database" else {}))
        user_mol = source in ("single", "custom")
        score_lo = 45.0 if user_mol else 32.0

        # 1. prepare complex
        self.progress("prepare", "Preparing the protein-ligand complex", 3)
        self.run_step([py, self.py("prepare_complex.py"), receptor, selector, "-o", "PROT.pdb",
                       "--json", "prepare.json"],
                      "The structure could not be prepared.", EXIT_LIGAND_NOT_FOUND,
                      code_map={2: EXIT_INPUT_MISSING, 3: EXIT_PROTEIN_SHORT, 13: EXIT_LIGAND_NOT_FOUND})
        info = read_json(self.p("prepare.json"), {}) or {}
        resname = (info.get("ligand") or {}).get("resname") or selector.split(":")[0].split()[0]

        # 2. tleap
        self.progress("tleap", "Standardising the structure with AmberTools tleap", 6)
        tleap = find_tleap()
        env = os.environ.copy()
        if tleap and not env.get("AMBERHOME"):
            env["AMBERHOME"] = os.path.dirname(os.path.dirname(tleap))
        shutil.copyfile(self.py("leap.cmd"), self.p("leap.cmd"))
        rm_f(self.p("PROT_LIG.pdb"))
        out = self.run_step([tleap or "tleap", "-f", "leap.cmd"],
                            "AmberTools tleap could not process the structure.", EXIT_TLEAP, timeout=1800,
                            hint="Check for unusual residues, missing backbone atoms or alternate locations.",
                            env=env)
        if not os.path.isfile(self.p("PROT_LIG.pdb")) or os.path.getsize(self.p("PROT_LIG.pdb")) == 0:
            raise PipelineError("AmberTools tleap did not write PROT_LIG.pdb.", EXIT_TLEAP,
                                details=tail_text(self.p("leap.log"), 30) or out)

        # 3. protein + pocket features
        self.progress("features", "Computing protein and pocket descriptors", 15)
        self.run_step([py, self.py("protein_features.py"), "PROT_LIG.pdb", resname],
                      "Protein / pocket descriptors could not be calculated.", EXIT_PROTEIN_FEATURES)
        self.progress("pocket", "Estimating pocket volume and surface area", 22)
        self.run_step([py, self.py("pocket_features.py"), "PROT_LIG.pdb", resname],
                      "Pocket volume / SASA could not be calculated.", EXIT_POCKET_FEATURES)
        self.run_step([py, self.py("merge_protein_features.py")],
                      "Protein descriptors could not be merged.", EXIT_MERGE)

        # 4. user molecules
        if source == "single":
            self.progress("ligand", "Building a 3D model of your molecule", 33)
            self.run_step([py, self.py("convert_mol.py"), os.path.basename(mol_file), "Ligand.sdf"],
                          "Your molecule could not be read or converted to 3D.", EXIT_MOLECULE,
                          hint="Check the SMILES / file format; only one molecule is allowed.")
            self.progress("ligand_features", "Computing ligand descriptors", 40)
            self.run_step([py, self.py("ligand_features.py"), "-i", "Ligand.sdf", "-o", "dataset.csv"],
                          "Descriptors could not be calculated for your molecule.", EXIT_NO_VALID_MOLECULES)
            mol_inputs = [self.p("dataset.csv")]
        elif source == "custom":
            self.progress("ligand_features", "Preparing and featurising your molecules", 33)
            max_mols = int((self.spec.get("options") or {}).get("max_molecules") or 0)
            cmd = [py, self.py("ligand_features.py"), "-i", os.path.basename(mol_file), "-o", "dataset.csv"]
            if max_mols:
                cmd += ["--max-molecules", str(max_mols)]
            self.run_step(cmd, "None of the uploaded molecules could be processed.", EXIT_NO_VALID_MOLECULES,
                          hint="Upload a valid multi-molecule SDF file.")
            mol_inputs = [self.p("dataset.csv")]

        # 5. scoring
        return self._score(mol_inputs, self.p("protein.csv"), "ligand", params, known_rows,
                           (score_lo, 98.0), large=(db_entry or {}).get("size") == "large")

    # -- reverse screening -------------------------------------------------
    def _reverse(self) -> dict:
        r = self.spec.get("reverse") or {}
        query = self._input(r.get("query_file"), EXIT_INPUT_MISSING, "query molecule")
        lib = resolve_catalog_key("libraries", r.get("library") or "homo")
        if not lib:
            raise PipelineError(f"Unknown target library '{r.get('library')}'.", EXIT_BAD_MODE)
        inputs, known = resolve_dataset_inputs(lib, "libraries")
        py = sys.executable
        self.progress("ligand", "Building a 3D model of your molecule", 3)
        self.run_step([py, self.py("convert_mol.py"), os.path.basename(query), "Ligand.sdf"],
                      "Your molecule could not be read or converted to 3D.", EXIT_MOLECULE,
                      hint="Check the SMILES / file format; only one molecule is allowed.")
        self.progress("ligand_features", "Computing ligand descriptors", 12)
        self.run_step([py, self.py("ligand_features.py"), "-i", "Ligand.sdf", "-o", "dataset.csv"],
                      "Descriptors could not be calculated for your molecule.", EXIT_NO_VALID_MOLECULES)
        return self._score(inputs, self.p("dataset.csv"), "protein", None, known, (20.0, 98.0),
                           large=lib.get("size") == "large")

    # -- scoring step ------------------------------------------------------
    def _score(self, inputs: List[str], fixed_csv: str, role: str, params: Optional[str],
               known_rows: Optional[int], prange: Tuple[float, float], large: bool) -> dict:
        import scoring

        opts = self.spec.get("options") or {}
        top = int(opts.get("top") or 0)
        try:
            total = known_rows or scoring.estimate_total_rows(scoring.expand_inputs(inputs))[0]
        except scoring.ScoringError as exc:
            raise PipelineError(str(exc), exc.exit_code)
        workers = self.workers
        threads = self.threads
        if not workers:
            workers = self.cpus if (total or 0) > 200_000 and self.cpus > 1 else 1
        if workers == 1 and not threads:
            threads = self.cpus
        label = "Scoring molecules" if role == "ligand" else "Scoring protein targets"
        self.progress("scoring", label, prange[0], f"Loading the ensemble models ({total:,} rows to score)")
        self.log(f"Scoring with {workers} worker(s), {threads or 1} thread(s) each; top limit: {top or 'none'}")
        try:
            return scoring.run_screening(
                inputs=inputs, fixed_csv=fixed_csv, models_dir=models_dir(), output=self.p("results.csv"),
                role=role, params_file=params, preview=self.p("results_preview.csv"),
                preview_rows=int(opts.get("preview_rows") or 2000), summary=self.p("summary.json"),
                progress=self.p("progress.json"), progress_range=prange, workers=workers, threads=threads,
                top=top, total_rows=known_rows, chunksize=int(opts.get("chunksize") or 50_000),
            )
        except scoring.ScoringError as exc:
            raise PipelineError(str(exc), exc.exit_code)

    # -- finalize ----------------------------------------------------------
    def _finalize(self, summary: dict) -> None:
        self.progress("finalize", "Finalising results", 99)
        res = self.p("results.csv")
        if not os.path.isfile(res):
            raise PipelineError("Scoring finished but results.csv was not written.", EXIT_INTERNAL)
        threshold = int(float(os.environ.get("SEARCHML_GZIP_MB", "50")) * 1024 * 1024)
        if os.path.getsize(res) > threshold:
            self.log("Compressing results.csv -> results.csv.gz")
            with open(res, "rb") as src, gzip.open(self.p("results.csv.gz.tmp"), "wb", compresslevel=6) as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
            os.replace(self.p("results.csv.gz.tmp"), self.p("results.csv.gz"))
        if summary.get("rows_written", 0) == 0:
            self.log("Note: no rows passed the filters - results.csv contains only the header.")

    # -- simulation (plumbing tests without RDKit / AmberTools / models) ---
    def _simulate(self, mode: str) -> dict:
        import scoring  # constants only

        opts = self.spec.get("options") or {}
        fail_at = opts.get("simulate_fail")
        speed = float(opts.get("simulate_seconds") or 6)
        stages = ([("prepare", "Preparing the protein-ligand complex", 3),
                   ("tleap", "Standardising the structure with AmberTools tleap", 6),
                   ("features", "Computing protein and pocket descriptors", 15),
                   ("pocket", "Estimating pocket volume and surface area", 22)]
                  if mode == "forward" else
                  [("ligand", "Building a 3D model of your molecule", 3),
                   ("ligand_features", "Computing ligand descriptors", 12)])
        for stage, label, pct in stages:
            self.progress(stage, label, pct)
            time.sleep(speed * 0.08)
            if self.cancel_requested:
                raise Cancelled()
            if fail_at == stage:
                raise PipelineError("Simulated failure at stage '%s'." % stage, EXIT_INTERNAL,
                                    details="(simulate_fail option set in job.json)",
                                    hint="This is a test job.")
        role = "ligand" if mode == "forward" else "protein"
        ids = self._simulated_ids(mode)
        total = len(ids)
        rng = random.Random(os.path.basename(self.dir))
        lo = 32.0 if mode == "forward" else 20.0
        rows = []
        prog = scoring.Progress(self.p("progress.json"), total, False,
                                "Scoring molecules" if role == "ligand" else "Scoring protein targets",
                                (lo, 98.0), interval=0.2)
        steps = 20
        for s in range(steps):
            time.sleep(speed * 0.6 / steps)
            if self.cancel_requested:
                raise Cancelled()
            done = int(total * (s + 1) / steps)
            prog.update(done, done, force=True)
        stats = scoring.SummaryStats()
        import pandas as pd
        for ident in ids:
            pki = max(2.0, min(11.5, rng.gauss(5.6, 1.1)))
            base = {k: pki + rng.gauss(0, 0.35) for k in ("xgb", "cat", "rf", "lgbm")}
            rows.append({"_pki": pki, "PDB_File": ident,
                         scoring.META_COLS[0]: round(pki, 2),
                         scoring.META_COLS[1]: round(scoring.KCAL_PER_PKI * pki, 2),
                         scoring.META_COLS[2]: round(10 ** (9 - pki), 2),
                         **{scoring.MODEL_COLS[k]: round(scoring.KCAL_PER_PKI * v, 2) for k, v in base.items()}})
        df = pd.DataFrame(rows, columns=["_pki", "PDB_File"] + scoring.RESULT_COLS)
        stats.update(df, "PDB_File")
        df = df.sort_values("_pki", ascending=False, kind="mergesort")
        top = int(opts.get("top") or 0)
        if top:
            df = df.head(top)
        out = df.drop(columns=["_pki"])
        out.to_csv(self.p("results.csv"), index=False)
        out.head(int(opts.get("preview_rows") or 2000)).to_csv(self.p("results_preview.csv"), index=False)
        summary = {"engine_version": scoring.__version__, "role": role, "id_column": "PDB_File",
                   "rows_read": total, "filtered_out": 0, "invalid_rows": 0, "rows_written": len(out),
                   "top_limit": top or None, "filters": {}, "model_order": scoring.DEFAULT_MODEL_ORDER,
                   "simulated": True, "workers": 1}
        summary.update(stats.as_dict())
        atomic_write_json(self.p("summary.json"), summary)
        return summary

    def _simulated_ids(self, mode: str) -> List[str]:
        n = int((self.spec.get("options") or {}).get("simulate_rows") or 0)
        try:
            if mode == "forward":
                f = self.spec.get("forward") or {}
                if (f.get("source") or "database") == "database":
                    entry = resolve_catalog_key("databases", f.get("database", ""))
                    if entry:
                        files, _ = resolve_dataset_inputs(entry, "databases")
                        import scoring
                        files = scoring.expand_inputs(files)
                        import pandas as pd
                        ids = pd.read_csv(files[0], usecols=[0], nrows=n or 5000).iloc[:, 0].astype(str)
                        return ids.tolist()
                elif f.get("source") == "single":
                    return ["query_molecule"]
                else:
                    count = 0
                    try:
                        with open(self.p(f.get("molecule_file", "")), "r", errors="replace") as fh:
                            count = sum(1 for ln in fh if ln.startswith("$$$$"))
                    except OSError:
                        pass
                    return [f"Molecule_{i + 1}" for i in range(max(1, count))]
            else:
                lib = resolve_catalog_key("libraries", (self.spec.get("reverse") or {}).get("library") or "homo")
                if lib:
                    files, _ = resolve_dataset_inputs(lib, "libraries")
                    import scoring
                    import pandas as pd
                    files = scoring.expand_inputs(files)
                    return pd.read_csv(files[0], usecols=[0], nrows=n or 5000).iloc[:, 0].astype(str).tolist()
        except (PipelineError, OSError, ValueError):
            pass
        prefix = "SIM" if mode == "forward" else "1ABC_LIG_A_"
        return [f"{prefix}{i:05d}" for i in range(1, (n or 500) + 1)]


_ERR_RE = re.compile(r"^\s*(?:\[!\]\s*)?(?:error|fatal)\s*:\s*(.+)$", re.IGNORECASE)


def _last_error_line(lines) -> str:
    for line in reversed(list(lines)):
        m = _ERR_RE.match(line.strip())
        if m:
            text = m.group(1).strip()
            return text if text.endswith(".") else text + "."
    return ""


def _shell_quote(s: str) -> str:
    return s if re.fullmatch(r"[A-Za-z0-9_./:=+-]+", s) else "'" + s.replace("'", "'\\''") + "'"


# --------------------------------------------------------------------------
# Command-line front ends (used by master_screening*.sh)
# --------------------------------------------------------------------------

MODE_ALIASES = {"DB": "drugbank", "DRUGBANK": "drugbank", "FDA": "fda", "BIMP": "bimp"}


def _new_job_dir(base: Optional[str]) -> str:
    if base:
        os.makedirs(base, exist_ok=False)
        return os.path.realpath(base)
    stamp = str(int(time.time()))
    path = os.path.join(os.getcwd(), stamp)
    n = 1
    while os.path.exists(path):
        path = os.path.join(os.getcwd(), f"{stamp}_{n}")
        n += 1
    os.makedirs(path)
    return path


def _autodetect(base: str, exts: List[str]) -> Optional[str]:
    for ext in exts:
        if os.path.isfile(base + ext):
            return base + ext
    return None


def cmd_forward(a) -> int:
    if not os.path.isfile(a.pdb):
        print(f"Error: PDB file '{a.pdb}' not found.")
        return EXIT_INPUT_MISSING
    mode = a.mode.upper()
    spec_f: Dict[str, object] = {"receptor": "receptor.pdb", "ligand": a.ligand}
    mol_src = None
    if mode == "SINGLE":
        mol_src = a.molecules or _autodetect("single_molecule", [".sdf", ".pdb", ".smi", ".mol", ".mol2"])
        if not mol_src or not os.path.isfile(mol_src):
            print("Error: 'SINGLE' mode selected, but no valid molecule file found.")
            print("Provide --molecules FILE or place single_molecule.{sdf,pdb,smi,mol} in this directory.")
            return EXIT_SINGLE_MISSING
        spec_f.update(source="single", molecule_file="single_molecule" + os.path.splitext(mol_src)[1].lower())
    elif mode == "CUSTOM":
        mol_src = a.molecules or ("custom_molecules.sdf" if os.path.isfile("custom_molecules.sdf") else None)
        if not mol_src or not os.path.isfile(mol_src):
            print("Error: 'CUSTOM' mode selected, but 'custom_molecules.sdf' not found.")
            return EXIT_CUSTOM_MISSING
        if not mol_src.lower().endswith(".sdf"):
            print("Error: custom mode supports .sdf files only.")
            return EXIT_CUSTOM_MISSING
        spec_f.update(source="custom", molecule_file="custom_molecules.sdf")
    else:
        key = MODE_ALIASES.get(mode, a.mode)
        entry = resolve_catalog_key("databases", key)
        if not entry:
            names = ", ".join(sorted({*MODE_ALIASES, *(e["key"].upper() for e in load_catalog("databases").values())}))
            print(f"Error: invalid mode '{a.mode}'. Allowed: {names}, SINGLE, CUSTOM")
            return EXIT_BAD_MODE
        spec_f.update(source="database", database=entry["key"])

    if a.no_filters:
        spec_f["filters"] = {}
    elif a.params:
        from scoring import parse_parameter_file
        spec_f["filters"] = {k: list(v) for k, v in parse_parameter_file(a.params).items() if k in FILTERABLE}

    job = _new_job_dir(a.job_dir)
    shutil.copyfile(a.pdb, os.path.join(job, "receptor.pdb"))
    if mol_src:
        shutil.copyfile(mol_src, os.path.join(job, str(spec_f["molecule_file"])))
    with open(os.path.join(job, "args.txt"), "w") as fh:
        fh.write(f"{a.pdb} , {a.ligand} , {a.mode}\n")
    spec = {"version": 2, "mode": "forward", "created_at": now_iso(), "forward": spec_f,
            "options": {"top": a.top}}
    atomic_write_json(os.path.join(job, "job.json"), spec)
    print(f"Job directory: {job}")
    return JobRunner(job, workers=a.workers, threads=a.threads, echo=True).execute()


def cmd_reverse(a) -> int:
    if a.smiles:
        src = None
    else:
        src = a.query or _autodetect("query_ligand", [".sdf", ".pdb", ".smi", ".mol", ".mol2"])
        if not src or not os.path.isfile(src):
            print("Error: no query molecule found. Expected query_ligand.sdf, .pdb or .smi "
                  "(or use --query FILE / --smiles STRING).")
            return EXIT_INPUT_MISSING
    lib = resolve_catalog_key("libraries", a.library)
    if not lib:
        names = ", ".join(sorted(e["key"] for e in load_catalog("libraries").values()))
        print(f"Error: unknown target library '{a.library}'. Available: {names}")
        return EXIT_BAD_MODE
    job = _new_job_dir(a.job_dir)
    if a.smiles:
        qname = "query_ligand.smi"
        with open(os.path.join(job, qname), "w") as fh:
            fh.write(a.smiles.strip() + " query_ligand\n")
    else:
        qname = "query_ligand" + os.path.splitext(src)[1].lower()
        shutil.copyfile(src, os.path.join(job, qname))
    spec = {"version": 2, "mode": "reverse", "created_at": now_iso(),
            "reverse": {"query_file": qname, "library": lib["key"]}, "options": {"top": a.top}}
    atomic_write_json(os.path.join(job, "job.json"), spec)
    print(f"Job directory: {job}")
    return JobRunner(job, workers=a.workers, threads=a.threads, echo=True).execute()


def cmd_run(a) -> int:
    try:
        runner = JobRunner(a.job_dir, workers=a.workers, threads=a.threads, echo=not a.quiet,
                           simulate=a.simulate)
    except PipelineError as exc:
        print(f"Error: {exc}")
        return exc.exit_code
    return runner.execute()


def cmd_check(a) -> int:
    problems = check_environment("forward", verbose=True)
    if problems:
        print("\nProblems found:")
        for p in problems:
            print("  - " + p)
        return EXIT_ENVIRONMENT
    print("\nEnvironment OK.")
    return 0


def cmd_catalog(a) -> int:
    out = {}
    for kind in ("databases", "libraries"):
        out[kind] = []
        for e in load_catalog(kind).values():
            try:
                files, count = resolve_dataset_inputs(e, kind)
                available = True
            except PipelineError:
                files, count, available = [], None, False
            out[kind].append({"key": e["key"], "label": e["label"], "available": available,
                              "rows": count, "size": e["size"], "path": files[0] if files else e["path"]})
    if a.json:
        print(json.dumps(out, indent=1))
    else:
        for kind, rows in out.items():
            print(f"\n{kind.capitalize()}:")
            for r in rows:
                state = "installed" if r["available"] else "missing"
                print(f"  {r['key']:<16} {r['label']:<40} {state}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="SEARCH-ML pipeline orchestrator",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--version", action="version", version=f"SEARCH-ML pipeline {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--workers", type=int, default=None, help="Parallel scoring processes")
        p.add_argument("--threads", type=int, default=None, help="Threads per model")
        p.add_argument("--job-dir", default=None, help="Job directory to create (default: timestamp)")
        p.add_argument("--top", type=int, default=0, help="Keep only the N best hits (0 = all)")

    pf = sub.add_parser("forward", help="Protein -> molecule library screening")
    pf.add_argument("--pdb", required=True, help="Protein-ligand complex (.pdb)")
    pf.add_argument("--ligand", required=True, help="Reference ligand: RES, RES:CHAIN:RESNUM or 'RES C NNNN'")
    pf.add_argument("--mode", required=True, help="DB, FDA, BIMP, any library key, SINGLE or CUSTOM")
    pf.add_argument("--molecules", help="Molecule file for SINGLE / CUSTOM modes")
    pf.add_argument("--params", help="Property filter file (Name: min,max)")
    pf.add_argument("--no-filters", action="store_true", help="Do not apply property filters")
    common(pf)
    pf.set_defaults(func=cmd_forward)

    pr = sub.add_parser("reverse", help="Molecule -> protein target library screening")
    pr.add_argument("--query", help="Query molecule (.sdf/.mol/.pdb/.smi)")
    pr.add_argument("--smiles", help="Query molecule as SMILES")
    pr.add_argument("--library", default="homo", help="Target library key (default: homo)")
    common(pr)
    pr.set_defaults(func=cmd_reverse)

    prun = sub.add_parser("run", help="Execute a prepared job directory")
    prun.add_argument("job_dir")
    prun.add_argument("--workers", type=int, default=None)
    prun.add_argument("--threads", type=int, default=None)
    prun.add_argument("--simulate", action="store_true", help="Test mode: fake results, no chemistry")
    prun.add_argument("--quiet", action="store_true", help="Do not echo the log to stdout")
    prun.set_defaults(func=cmd_run)

    pc = sub.add_parser("check", help="Verify the environment")
    pc.set_defaults(func=cmd_check)

    pcat = sub.add_parser("catalog", help="List compound and target libraries")
    pcat.add_argument("--json", action="store_true")
    pcat.set_defaults(func=cmd_catalog)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
