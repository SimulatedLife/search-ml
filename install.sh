#!/usr/bin/env bash
# ==============================================================================
# SEARCH-ML : installer for a workstation or a single server
#
# Creates (or updates) the "search-ml" conda environment, records SEARCH_ML_HOME
# in it, fetches the Random Forest model if it is missing, and checks that
# everything works.
#
# Usage:
#   ./install.sh                     interactive install
#   ./install.sh --yes               no questions (download the model if missing)
#   ./install.sh --name myenv        use another environment name
#   ./install.sh --skip-env          only check / download models
#
# For a cluster that runs jobs submitted from the SEARCH-ML web portal use
# hpc/install_hpc.sh instead (it calls this script).
#
# Author : Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi (dheeraj@scfbio-iitd.res.in)
# ==============================================================================
set -Eeuo pipefail

ENV_NAME="search-ml"
ASSUME_YES=0
SKIP_ENV=0
RF_URL="${SEARCHML_RF_URL:-https://scfbio.iitd.ac.in/search-ml/models/final_rf_model.joblib}"
RF_SHA256="8d8573fd152739df1b2d25e174c63e1314412dd6df3a7b61a56bab1c571fe8e0"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [ $# -gt 0 ]; do
  case "$1" in
    -y|--yes) ASSUME_YES=1 ;;
    -n|--name) ENV_NAME="$2"; shift ;;
    --skip-env) SKIP_ENV=1 ;;
    -h|--help) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option: $1 (use --help)"; exit 2 ;;
  esac
  shift
done

say()  { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
ok()   { printf '\033[1;32m ok\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m !!\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31mERR\033[0m %s\n' "$*" >&2; exit 1; }
ask()  { [ "$ASSUME_YES" -eq 1 ] && return 0; read -r -p "$1 [Y/n] " r; [[ -z "$r" || "$r" =~ ^[Yy] ]]; }
trap 'die "Installation failed at line $LINENO. Fix the problem above and re-run ./install.sh"' ERR

cd "$ROOT"
[ -f environment.yml ] || die "environment.yml not found - run this script from the search-ml folder."

# ------------------------------------------------------------------ conda ----
detect_conda() {
  local c
  for c in "${CONDA_EXE:-}" "$(command -v conda 2>/dev/null || true)" \
           "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" "$HOME/miniconda3/bin/conda" \
           "$HOME/anaconda3/bin/conda" /opt/conda/bin/conda /opt/miniconda3/bin/conda \
           /opt/miniforge3/bin/conda /usr/local/miniconda3/bin/conda; do
    if [ -n "$c" ] && [ -x "$c" ]; then
      CONDA_BASE="$("$c" info --base 2>/dev/null)" && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ] && return 0
    fi
  done
  return 1
}

if [ "$SKIP_ENV" -eq 0 ]; then
  say "Looking for conda"
  if ! detect_conda; then
    cat <<MSG
Conda (Miniforge, Miniconda or Anaconda) is required but was not found.
Install Miniforge (recommended, uses conda-forge by default):
   curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-\$(uname)-\$(uname -m).sh
   bash Miniforge3-\$(uname)-\$(uname -m).sh
then open a new terminal and run ./install.sh again.
MSG
    exit 1
  fi
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
  ok "conda at $CONDA_BASE"
  SOLVER="conda"
  command -v mamba >/dev/null 2>&1 && SOLVER="mamba"

  if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    say "Updating environment '$ENV_NAME' (this can take a while)"
    $SOLVER env update -n "$ENV_NAME" -f environment.yml --prune
  else
    say "Creating environment '$ENV_NAME' (AmberTools is large: 10-30 min)"
    $SOLVER env create -n "$ENV_NAME" -f environment.yml
  fi
  conda env config vars set SEARCH_ML_HOME="$ROOT" -n "$ENV_NAME" >/dev/null
  ok "SEARCH_ML_HOME=$ROOT recorded in the environment"
  conda activate "$ENV_NAME"
fi

# ------------------------------------------------------------------ models ---
say "Checking model files"
RF="$ROOT/models/final_rf_model.joblib"
if [ ! -s "$RF" ]; then
  warn "The Random Forest model (~830 MB) is not in models/ (too large for GitHub)."
  if ask "Download it now from $RF_URL ?"; then
    if command -v curl >/dev/null 2>&1; then
      curl -fL --retry 3 -C - -o "$RF.part" "$RF_URL"
    else
      wget -c -O "$RF.part" "$RF_URL"
    fi
    mv "$RF.part" "$RF"
  else
    warn "Download it later from https://scfbio.iitd.ac.in/search-ml/download.php into models/"
  fi
fi
if [ -s "$RF" ] && command -v sha256sum >/dev/null 2>&1; then
  if echo "$RF_SHA256  $RF" | sha256sum -c --status; then ok "final_rf_model.joblib checksum verified"
  else warn "final_rf_model.joblib checksum does not match - the download may be incomplete."; fi
fi
if command -v sha256sum >/dev/null 2>&1; then ( cd models && sha256sum -c --quiet --ignore-missing SHA256SUMS ) && ok "model checksums match"; fi

# ------------------------------------------------------------------- check ---
chmod +x master_screening.sh master_screening_reverse.sh scripts/*.py hpc/bin/* 2>/dev/null || true
if [ "$SKIP_ENV" -eq 0 ]; then
  say "Verifying the environment"
  if SEARCH_ML_HOME="$ROOT" python scripts/pipeline.py check; then
    ok "SEARCH-ML is ready"
  else
    warn "Some checks failed (see above)."
  fi
fi

cat <<MSG

Next steps
  conda activate $ENV_NAME
  cd example
  ../master_screening.sh 9kte.pdb 9UM FDA          # forward screen against FDA drugs
  ../master_screening_reverse.sh -s "CC(=O)Oc1ccccc1C(=O)O"   # target identification

If 'echo \$SEARCH_ML_HOME' prints nothing after activation, run:
  conda env config vars set SEARCH_ML_HOME="$ROOT" -n $ENV_NAME && conda deactivate && conda activate $ENV_NAME
MSG
