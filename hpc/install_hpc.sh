#!/usr/bin/env bash
# ==============================================================================
# SEARCH-ML : HPC installer
#
# Sets up a cluster (or a single Linux server) to run the jobs submitted from the
# SEARCH-ML web portal:
#   1. creates the conda environment and checks the models   (../install.sh)
#   2. detects the batch system (SLURM / PBS / none) and writes hpc/searchml-hpc.conf
#   3. creates the jobs directory
#   4. optionally authorises the web server's SSH key, restricted to SEARCH-ML commands
#   5. runs a self-test job and prints the settings to paste into the web portal config
#
# Usage:
#   bash hpc/install_hpc.sh                       # interactive
#   bash hpc/install_hpc.sh --yes --scheduler slurm --jobs-dir /scratch/searchml/jobs \
#        --partition compute --web-key "ssh-ed25519 AAAA... searchml-web"
#
# Options:
#   --yes                   accept defaults, no questions
#   --scheduler NAME        local | slurm | pbs   (default: auto-detect)
#   --jobs-dir DIR          job folders (must be visible from the compute nodes)
#   --partition NAME        SLURM partition / PBS queue
#   --account NAME          SLURM/PBS account
#   --env NAME              conda environment name (default: search-ml)
#   --web-key "KEY"         public key of the web server to authorise (forced command)
#   --skip-env              do not create/update the conda environment
#   --real-test             also run a real test job (needs the RF model; ~1-3 min)
#   --force                 overwrite an existing hpc/searchml-hpc.conf
#
# Author : Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi (dheeraj@scfbio-iitd.res.in)
# ==============================================================================
set -Eeuo pipefail

HPC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HPC_DIR")"
CONF="$HPC_DIR/searchml-hpc.conf"

YES=0 SCHED="" JOBS="" PART="" ACCT="" ENV_NAME="search-ml" WEBKEY="" SKIP_ENV=0 REAL_TEST=0 FORCE=0
while [ $# -gt 0 ]; do
  case "$1" in
    --yes|-y) YES=1 ;;
    --scheduler) SCHED="$2"; shift ;;
    --jobs-dir) JOBS="$2"; shift ;;
    --partition) PART="$2"; shift ;;
    --account) ACCT="$2"; shift ;;
    --env) ENV_NAME="$2"; shift ;;
    --web-key) WEBKEY="$2"; shift ;;
    --skip-env) SKIP_ENV=1 ;;
    --real-test) REAL_TEST=1 ;;
    --force) FORCE=1 ;;
    -h|--help) sed -n '2,32p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option $1 (see --help)"; exit 2 ;;
  esac
  shift
done

say()  { printf '\n\033[1;34m==>\033[0m \033[1m%s\033[0m\n' "$*"; }
ok()   { printf '\033[1;32m ok\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m !!\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31mERR\033[0m %s\n' "$*" >&2; exit 1; }
ask()  { # ask "Question" default -> echoes answer
  local q="$1" d="${2:-}" r
  if [ "$YES" -eq 1 ]; then echo "$d"; return; fi
  read -r -p "$q [${d}]: " r </dev/tty || true
  echo "${r:-$d}"
}
yesno() { [ "$YES" -eq 1 ] && return 0; local r; read -r -p "$1 [Y/n] " r </dev/tty || true; [[ -z "$r" || "$r" =~ ^[Yy] ]]; }
trap 'die "HPC installation failed at line $LINENO."' ERR

echo "SEARCH-ML HPC installer - repository: $ROOT"

# 1 -------------------------------------------------------------------- env ---
if [ "$SKIP_ENV" -eq 0 ]; then
  say "Step 1/5: conda environment and models"
  bash "$ROOT/install.sh" --yes --name "$ENV_NAME"
else
  say "Step 1/5: skipped (--skip-env)"
fi

CONDA_BASE=""
for c in "${CONDA_EXE:-}" "$(command -v conda 2>/dev/null || true)" "$HOME/miniforge3/bin/conda" \
         "$HOME/mambaforge/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
  if [ -n "$c" ] && [ -x "$c" ]; then CONDA_BASE="$("$c" info --base 2>/dev/null || true)"; [ -n "$CONDA_BASE" ] && break; fi
done
ENV_PY=""
[ -n "$CONDA_BASE" ] && [ -x "$CONDA_BASE/envs/$ENV_NAME/bin/python" ] && ENV_PY="$CONDA_BASE/envs/$ENV_NAME/bin/python"
if [ -z "$ENV_PY" ]; then
  warn "Could not find the '$ENV_NAME' environment python - set PYTHON= in $CONF later."
else
  ok "environment python: $ENV_PY"
fi

# 2 -------------------------------------------------------------- scheduler ---
say "Step 2/5: batch system"
if [ -z "$SCHED" ]; then
  if command -v sbatch >/dev/null 2>&1; then SCHED=slurm
  elif command -v qsub >/dev/null 2>&1; then SCHED=pbs
  else SCHED=local; fi
  SCHED="$(ask "Scheduler to use (local / slurm / pbs)" "$SCHED")"
fi
case "$SCHED" in local|slurm|pbs) ok "scheduler: $SCHED" ;; *) die "unknown scheduler '$SCHED'";; esac

if [ "$SCHED" = "slurm" ] && [ -z "$PART" ] && command -v sinfo >/dev/null 2>&1; then
  DEF_PART="$(sinfo -h -o '%P' 2>/dev/null | grep '\*' | tr -d '*' | head -n1 || true)"
  echo "    partitions: $(sinfo -h -o '%P' 2>/dev/null | tr '\n' ' ')"
  PART="$(ask "SLURM partition for SEARCH-ML jobs (empty = default)" "$DEF_PART")"
fi

# 3 ------------------------------------------------------------------- jobs ---
say "Step 3/5: jobs directory"
if [ -z "$JOBS" ]; then
  JOBS="$(ask "Directory for job folders (shared with compute nodes)" "$HOME/searchml-jobs")"
fi
JOBS="${JOBS/#\~/$HOME}"
mkdir -p "$JOBS"
chmod 750 "$JOBS" || true
ok "jobs directory: $JOBS ($(df -h "$JOBS" | awk 'NR==2{print $4" free"}'))"
if [ "$SCHED" != "local" ]; then
  warn "Make sure $JOBS and $ROOT are mounted on every compute node."
fi

if [ -f "$CONF" ] && [ "$FORCE" -eq 0 ]; then
  cp "$CONF" "$CONF.bak.$(date +%Y%m%d%H%M%S)"
  warn "Existing $CONF backed up; writing a new one."
fi
CPUS_LARGE=16
[ "$SCHED" = "local" ] && CPUS_LARGE="$(( $(nproc 2>/dev/null || echo 4) > 16 ? 16 : $(nproc 2>/dev/null || echo 4) ))"
sed -e "s|^SEARCH_ML_HOME=.*|SEARCH_ML_HOME=$ROOT|" \
    -e "s|^JOBS_DIR=.*|JOBS_DIR=$JOBS|" \
    -e "s|^CONDA_BASE=.*|CONDA_BASE=$CONDA_BASE|" \
    -e "s|^CONDA_ENV=.*|CONDA_ENV=$ENV_NAME|" \
    -e "s|^SCHEDULER=.*|SCHEDULER=$SCHED|" \
    -e "s|^PARTITION=.*|PARTITION=$PART|" \
    -e "s|^ACCOUNT=.*|ACCOUNT=$ACCT|" \
    -e "s|^LARGE_CPUS=.*|LARGE_CPUS=$CPUS_LARGE|" \
    "$HPC_DIR/searchml-hpc.conf.example" > "$CONF"
chmod 640 "$CONF"
ok "wrote $CONF (edit it to change CPUs, memory, time limits)"
chmod +x "$HPC_DIR"/bin/* "$ROOT"/master_screening*.sh "$ROOT"/scripts/*.py 2>/dev/null || true

# row counts for the progress bars
if [ -n "$ENV_PY" ]; then
  "$ENV_PY" "$ROOT/scripts/dataset_tools.py" count --catalog >/dev/null 2>&1 || true
fi

# 4 --------------------------------------------------------------- ssh key ---
say "Step 4/5: web server access"
GATEWAY="$HPC_DIR/bin/searchml-gateway"
if [ -z "$WEBKEY" ] && [ "$YES" -eq 0 ]; then
  echo "Paste the web server's SSH public key (one line, e.g. 'ssh-ed25519 AAAA... searchml-web'),"
  echo "or press Enter to skip and do it later:"
  read -r WEBKEY </dev/tty || true
fi
if [ -n "$WEBKEY" ]; then
  if ! [[ "$WEBKEY" =~ ^(ssh-(ed25519|rsa)|ecdsa-sha2-nistp(256|384|521))[[:space:]]+[A-Za-z0-9+/=]+ ]]; then
    die "That does not look like an SSH public key."
  fi
  mkdir -p "$HOME/.ssh" && chmod 700 "$HOME/.ssh"
  touch "$HOME/.ssh/authorized_keys" && chmod 600 "$HOME/.ssh/authorized_keys"
  KEYBODY="$(echo "$WEBKEY" | awk '{print $2}')"
  if grep -qF "$KEYBODY" "$HOME/.ssh/authorized_keys"; then
    warn "This key is already in ~/.ssh/authorized_keys - left unchanged."
  else
    echo "command=\"$GATEWAY\",no-port-forwarding,no-X11-forwarding,no-agent-forwarding,no-pty $WEBKEY" \
      >> "$HOME/.ssh/authorized_keys"
    ok "key authorised (restricted to SEARCH-ML commands via $GATEWAY)"
  fi
else
  echo "    Later, add this line to ~/.ssh/authorized_keys (replace the key):"
  echo "    command=\"$GATEWAY\",no-port-forwarding,no-X11-forwarding,no-agent-forwarding,no-pty ssh-ed25519 AAAA...  searchml-web"
fi

# 5 ---------------------------------------------------------------- tests ---
say "Step 5/5: self-test"
HPC_CMD="$HPC_DIR/bin/searchml-hpc"
if "$HPC_CMD" ping > /tmp/searchml_ping.$$ 2>&1; then
  if command -v python3 >/dev/null 2>&1; then
    python3 - "/tmp/searchml_ping.$$" <<'PY' || true
import json, sys
d = json.load(open(sys.argv[1]))
print("    host %s | scheduler %s | %.1f GB free | ready: %s" % (d["host"], d["scheduler"], d["disk_free_gb"], d["ready"]))
for p in d.get("problems", []):
    print("    problem: " + p)
inst = [k for k, v in d["databases"].items() if v["available"]]
miss = [k for k, v in d["databases"].items() if not v["available"]]
print("    compound libraries installed: %s" % (", ".join(inst) or "none"))
if miss:
    print("    not installed (greyed out on the web site): %s" % ", ".join(miss))
libs = [k for k, v in d["libraries"].items() if v["available"]]
print("    target libraries installed: %s" % (", ".join(libs) or "none"))
PY
  else
    cat /tmp/searchml_ping.$$
  fi
  ok "searchml-hpc ping works"
else
  cat /tmp/searchml_ping.$$; rm -f /tmp/searchml_ping.$$
  die "searchml-hpc ping failed"
fi
rm -f /tmp/searchml_ping.$$

if "$HPC_CMD" selftest --simulate >/dev/null 2>&1; then ok "simulated job ran end to end"
else warn "simulated self-test failed - run: $HPC_CMD selftest --simulate"; fi

if [ "$REAL_TEST" -eq 1 ] || { [ "$YES" -eq 0 ] && yesno "Run a real test job now (9KTE vs FDA drugs, needs the RF model, 1-3 min)?"; }; then
  if "$HPC_CMD" selftest --real > /tmp/searchml_selftest.$$ 2>&1; then
    ok "real test job completed:"; python3 -c "import json,sys;[print('    '+r) for r in json.load(open(sys.argv[1]))['top_rows']]" /tmp/searchml_selftest.$$ || true
  else
    warn "real test job failed:"; cat /tmp/searchml_selftest.$$
  fi
  rm -f /tmp/searchml_selftest.$$
fi

cat <<MSG

------------------------------------------------------------------------------
 SEARCH-ML HPC side is installed.

 Web portal settings (config/config.php on the web server):
     'hpc' => [
         'transport' => 'ssh2',            // or 'openssh'
         'host'      => '$(hostname -f 2>/dev/null || hostname)',
         'port'      => 22,
         'user'      => '$(whoami)',
         'key_file'  => '/etc/searchml/id_ed25519',   // private key readable by the web server
         'command'   => '$HPC_CMD',
         'jobs_dir'  => '$JOBS',
     ],

 Optional cron job on the HPC to delete finished jobs after the retention period:
     17 3 * * *  $HPC_CMD cleanup >/dev/null 2>&1

 Useful commands:
     $HPC_CMD ping            # health + installed libraries
     $HPC_CMD list            # recent jobs
     $HPC_CMD status <job>    # one job
     $HPC_CMD selftest --real # end-to-end test
------------------------------------------------------------------------------
MSG
