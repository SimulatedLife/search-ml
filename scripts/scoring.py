"""
SEARCH-ML scoring engine.

Shared by ``screening.py`` (one protein vs. many molecules) and
``reverse_screening.py`` (one molecule vs. many proteins).

What it does
------------
* Loads the stacked ensemble (scaler, XGBoost, CatBoost, Random Forest,
  LightGBM and the Ridge meta-model) once.
* Streams the large side of the screen (molecules or proteins) from one or
  more CSV files (``.csv`` or ``.csv.gz``) in chunks, so memory stays bounded
  even for 100M-row libraries.
* Applies the optional property filters (``parameter.txt``) *before*
  prediction, which gives exactly the same hits as filtering afterwards but
  skips work on molecules that would be discarded anyway.
* Optionally scores chunks in parallel worker processes (``workers > 1``).
* Writes a globally sorted ``results.csv`` (strongest predicted binders first)
  using an external merge sort, plus a small preview file and a
  ``summary.json`` with counts and an affinity histogram over *all* scored
  rows.
* Reports progress to a JSON file that the web portal / HPC runner polls.

Author: Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi (original scoring logic)
"""

from __future__ import annotations

import csv
import heapq
import json
import math
import multiprocessing as mp
import os
import re
import shutil
import sys
import tempfile
import time
from collections import deque
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

__version__ = "2.0.0"

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

#: pKi -> kcal/mol conversion used throughout SEARCH-ML (dG = -1.364 * pKi at 298 K)
KCAL_PER_PKI = -1.364

SCALER_FILE = "scaler_final.joblib"
META_FILE = "final_meta_model.joblib"
MODEL_FILES = {
    "xgb": "final_xgb_model.joblib",
    "cat": "final_cat_model.joblib",
    "rf": "final_rf_model.joblib",
    "lgbm": "final_lgbm_model.joblib",
}
ORDER_FILE = "final_ensemble_model_order.txt"

#: Order in which base-model predictions are stacked for the meta-model.
#: Must match parameters/final_ensemble_model_order.txt shipped with the models.
DEFAULT_MODEL_ORDER = ["cat", "xgb", "rf", "lgbm"]

META_COLS = ["PBE (pKi/pKd/pIC50)", "PBE (kcal/mol)", "PBE (nM)"]
MODEL_COLS = {
    "xgb": "XGB (kcal/mol)",
    "cat": "CAT (kcal/mol)",
    "rf": "RF (kcal/mol)",
    "lgbm": "LGBM (kcal/mol)",
}
RESULT_COLS = META_COLS + [MODEL_COLS[k] for k in ("xgb", "cat", "rf", "lgbm")]

LIGAND_ID_CANDIDATES = ["PDB_File", "Molecule", "Name", "ID", "ZINC_ID", "SMILES"]
PROTEIN_ID_CANDIDATES = ["Protein_ID", "Entry", "Accession", "Protein_Name", "PDB_ID", "PDB_File"]

#: Histogram of predicted pKi written to summary.json (edges 0.0 .. 15.0, step 0.1)
HIST_MIN, HIST_MAX, HIST_STEP = 0.0, 15.0, 0.1
HIST_BINS = int(round((HIST_MAX - HIST_MIN) / HIST_STEP))


class ScoringError(RuntimeError):
    """Raised for user-facing scoring problems (bad inputs, missing models...)."""

    def __init__(self, message: str, exit_code: int = 1):
        super().__init__(message)
        self.exit_code = exit_code


def log(msg: str) -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def clean_feature_name(name: str) -> str:
    """Make feature names compatible with XGBoost/LightGBM (same rule as training)."""
    name = re.sub(r"[,\(\)\-\s]", "_", str(name))
    name = re.sub(r"[\{\}\[\]:\"]", "", name)
    return name


def drop_unnamed(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in df.columns if "Unnamed" in str(c)]
    return df.drop(columns=cols) if cols else df


def select_id_column(columns: Sequence[str], candidates: Sequence[str]) -> str:
    for c in candidates:
        if c in columns:
            return c
    return columns[0]


def parse_parameter_file(file_path: Optional[str]) -> Dict[str, Tuple[float, float]]:
    """Parse ``Key: min,max`` (or ``Key: value``) lines. Missing file -> no filters."""
    params: Dict[str, Tuple[float, float]] = {}
    if not file_path:
        return params
    if not os.path.exists(file_path):
        log(f"Note: parameter file '{file_path}' not found - no property filters applied.")
        return params
    with open(file_path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            key, _, value = line.partition(":")
            key = key.strip()
            try:
                if "," in value:
                    lo, hi = value.split(",", 1)
                    lo_f, hi_f = float(lo.strip()), float(hi.strip())
                else:
                    lo_f = hi_f = float(value.strip())
            except ValueError:
                log(f"Warning: ignoring malformed filter line: {line!r}")
                continue
            if lo_f > hi_f:
                lo_f, hi_f = hi_f, lo_f
            params[key] = (lo_f, hi_f)
    return params


def find_order_file(models_dir: str, explicit: Optional[str] = None) -> Optional[str]:
    """Locate final_ensemble_model_order.txt (explicit path, models dir, cwd, $SEARCH_ML_HOME/parameters)."""
    candidates = []
    if explicit:
        candidates.append(explicit)
    candidates.append(os.path.join(models_dir, ORDER_FILE))
    candidates.append(os.path.join(os.getcwd(), ORDER_FILE))
    home = os.environ.get("SEARCH_ML_HOME")
    if home:
        candidates.append(os.path.join(home, "parameters", ORDER_FILE))
    here = os.path.dirname(os.path.realpath(__file__))
    candidates.append(os.path.join(os.path.dirname(here), "parameters", ORDER_FILE))
    for c in candidates:
        if c and os.path.isfile(c):
            return c
    return None


def read_model_order(path: Optional[str]) -> List[str]:
    if not path:
        log(f"Warning: {ORDER_FILE} not found - using default order {DEFAULT_MODEL_ORDER}.")
        return list(DEFAULT_MODEL_ORDER)
    with open(path, "r", encoding="utf-8") as fh:
        order = [ln.strip().lower() for ln in fh if ln.strip()]
    unknown = [m for m in order if m not in MODEL_FILES]
    if unknown or sorted(order) != sorted(MODEL_FILES):
        raise ScoringError(
            f"Invalid model order file {path}: expected the four keys {sorted(MODEL_FILES)}, got {order}"
        )
    return order


def atomic_write_json(path: str, payload: dict) -> None:
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, sort_keys=False, default=_json_default)
    os.replace(tmp, path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# --------------------------------------------------------------------------
# Model ensemble
# --------------------------------------------------------------------------

class Ensemble:
    """The stacked SEARCH-ML ensemble."""

    def __init__(self, models_dir: str, order_file: Optional[str] = None):
        import joblib  # imported lazily so --help works without the env

        self.models_dir = models_dir
        missing = [f for f in [SCALER_FILE, META_FILE, *MODEL_FILES.values()]
                   if not os.path.isfile(os.path.join(models_dir, f))]
        if missing:
            hint = ""
            if MODEL_FILES["rf"] in missing:
                hint = (" The Random Forest model (~830 MB) is not stored on GitHub; download it from "
                        "https://scfbio.iitd.ac.in/search-ml/download.php into the models/ folder.")
            raise ScoringError(f"Missing model file(s) in {models_dir}: {', '.join(missing)}.{hint}", 7)

        log(f"Loading models from '{models_dir}'...")
        t0 = time.time()
        self.scaler = joblib.load(os.path.join(models_dir, SCALER_FILE))
        self.models = {k: joblib.load(os.path.join(models_dir, f)) for k, f in MODEL_FILES.items()}
        self.meta = joblib.load(os.path.join(models_dir, META_FILE))
        self.order = read_model_order(find_order_file(models_dir, order_file))
        self.threads = None
        if hasattr(self.scaler, "feature_names_in_"):
            self.feature_columns = [str(c) for c in self.scaler.feature_names_in_]
        else:  # pragma: no cover - all shipped scalers carry feature names
            self.feature_columns = None
        log(f"All models loaded in {time.time() - t0:.1f} s (meta-model order: {', '.join(self.order)}).")

    # -- threading ---------------------------------------------------------
    def limit_threads(self, n: Optional[int]) -> None:
        """Limit the number of threads each model uses during prediction."""
        if not n or n < 1:
            return
        self.threads = int(n)
        rf = self.models.get("rf")
        if rf is not None and hasattr(rf, "n_jobs"):
            try:
                rf.n_jobs = self.threads
            except Exception:
                pass
        for key in ("xgb", "lgbm"):
            m = self.models.get(key)
            if m is not None and hasattr(m, "set_params"):
                try:
                    m.set_params(n_jobs=self.threads)
                except Exception:
                    pass

    # -- prediction --------------------------------------------------------
    def _predict_one(self, key: str, X: pd.DataFrame) -> np.ndarray:
        model = self.models[key]
        if self.threads:
            if key == "cat":
                try:
                    return np.asarray(model.predict(X, thread_count=self.threads), dtype=float)
                except TypeError:
                    pass
            elif key == "lgbm":
                try:
                    return np.asarray(model.predict(X, num_threads=self.threads), dtype=float)
                except TypeError:
                    pass
        return np.asarray(model.predict(X), dtype=float)

    def predict(self, X: pd.DataFrame) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """Return (meta pKi, {model: pKi}) for a feature frame in scaler column order."""
        scaled = pd.DataFrame(self.scaler.transform(X), columns=X.columns, index=X.index)
        base = {k: self._predict_one(k, scaled) for k in MODEL_FILES}
        stacked = np.column_stack([base[m] for m in self.order])
        meta = np.asarray(self.meta.predict(stacked), dtype=float)
        return meta, base


# --------------------------------------------------------------------------
# Input streaming
# --------------------------------------------------------------------------

def expand_inputs(paths: Sequence[str]) -> List[str]:
    """Expand directories into their CSV parts (sorted); keep files as given."""
    out: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            parts_dir = os.path.join(p, "parts") if os.path.isdir(os.path.join(p, "parts")) else p
            files = sorted(
                os.path.join(parts_dir, f) for f in os.listdir(parts_dir)
                if f.endswith((".csv", ".csv.gz")) and not f.startswith(".")
            )
            if not files:
                raise ScoringError(f"No .csv / .csv.gz files found in {parts_dir}", 2)
            out.extend(files)
        elif os.path.isfile(p):
            out.append(p)
        else:
            raise ScoringError(f"Input file not found: {p}", 2)
    return out


def _count_lines(path: str) -> int:
    n = 0
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(1 << 22)
            if not buf:
                break
            n += buf.count(b"\n")
    return n


def estimate_total_rows(paths: Sequence[str]) -> Tuple[int, bool]:
    """Return (rows, is_estimate). Uses <file>.count / count.txt caches when present."""
    total, estimated = 0, False
    seen_dirs = set()
    for p in paths:
        cache = p + ".count"
        d = os.path.dirname(p)
        dir_cache = os.path.join(d, "count.txt")
        if not os.path.isfile(dir_cache) and os.path.basename(d) == "parts":
            dir_cache = os.path.join(os.path.dirname(d), "count.txt")
        if os.path.isfile(cache):
            try:
                total += int(open(cache).read().split()[0])
                continue
            except (ValueError, IndexError):
                pass
        if os.path.isfile(dir_cache) and len(paths) > 1:
            if dir_cache not in seen_dirs:
                seen_dirs.add(dir_cache)
                try:
                    total += int(open(dir_cache).read().split()[0])
                except (ValueError, IndexError):
                    pass
            continue
        size = os.path.getsize(p)
        if p.endswith(".gz"):
            # rough: CSV feature rows compress ~4x, ~450 bytes per row uncompressed
            total += max(1, int(size * 4 / 450))
            estimated = True
        elif size <= 256 * 1024 * 1024:
            total += max(0, _count_lines(p) - 1)
        else:
            with open(p, "rb") as fh:
                sample = fh.read(1 << 20)
            lines = max(1, sample.count(b"\n"))
            total += int(size / (len(sample) / lines))
            estimated = True
    return total, estimated


def iter_chunks(paths: Sequence[str], chunksize: int) -> Iterator[pd.DataFrame]:
    for p in paths:
        try:
            reader = pd.read_csv(p, chunksize=chunksize, low_memory=False)
        except pd.errors.EmptyDataError:
            log(f"Warning: {p} is empty - skipped.")
            continue
        for chunk in reader:
            yield chunk


# --------------------------------------------------------------------------
# Chunk scoring (runs in the main process or in worker processes)
# --------------------------------------------------------------------------

_STATE: Dict[str, object] = {}


def _setup_state(ensemble: Ensemble, fixed: pd.DataFrame, filters: Dict[str, Tuple[float, float]],
                 id_col: str, nan_policy: str = "auto") -> None:
    feature_cols = ensemble.feature_columns
    fixed = drop_unnamed(fixed).copy()
    fixed.columns = [clean_feature_name(c) for c in fixed.columns]
    _STATE.clear()
    _STATE.update(
        ensemble=ensemble,
        fixed=fixed.iloc[0],
        filters=dict(filters),
        id_col=id_col,
        feature_cols=feature_cols,
        nan_policy=nan_policy,
    )


def _init_worker(threads: int, spawn_args: Optional[dict]) -> None:
    """Pool initializer. With 'fork' the state is inherited; with 'spawn' we rebuild it."""
    import signal as _signal
    _signal.signal(_signal.SIGTERM, _signal.SIG_DFL)  # let the parent handle cancellation
    _signal.signal(_signal.SIGINT, _signal.SIG_IGN)
    os.environ.setdefault("OMP_NUM_THREADS", str(threads))
    if spawn_args and "ensemble" not in _STATE:
        ens = Ensemble(spawn_args["models_dir"], spawn_args.get("order_file"))
        _setup_state(ens, spawn_args["fixed"], spawn_args["filters"], spawn_args["id_col"],
                     spawn_args.get("nan_policy", "auto"))
    ens = _STATE.get("ensemble")
    if ens is not None:
        ens.limit_threads(threads)


def score_chunk(chunk: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Filter, featurise and score one chunk. Returns (results, counters)."""
    ens: Ensemble = _STATE["ensemble"]  # type: ignore[assignment]
    fixed: pd.Series = _STATE["fixed"]  # type: ignore[assignment]
    filters: Dict[str, Tuple[float, float]] = _STATE["filters"]  # type: ignore[assignment]
    id_col: str = _STATE["id_col"]  # type: ignore[assignment]
    feature_cols: List[str] = _STATE["feature_cols"]  # type: ignore[assignment]

    counters = {"read": int(len(chunk)), "filtered": 0, "invalid": 0, "incomplete": 0, "scored": 0}
    chunk = drop_unnamed(chunk)
    if id_col not in chunk.columns:
        raise ScoringError(f"Identifier column '{id_col}' missing from input chunk.", 2)

    # 1. property filters on the raw descriptor columns
    if filters:
        mask = np.ones(len(chunk), dtype=bool)
        for name, (lo, hi) in filters.items():
            if name in chunk.columns:
                vals = pd.to_numeric(chunk[name], errors="coerce").to_numpy(dtype=float)
                with np.errstate(invalid="ignore"):
                    mask &= (vals >= lo) & (vals <= hi)
        counters["filtered"] = int((~mask).sum())
        if not mask.all():
            chunk = chunk.loc[mask]
    if chunk.empty:
        return _empty_result(id_col), counters

    # 2. assemble the feature matrix (chunk columns + broadcast fixed row)
    cleaned = {clean_feature_name(c): c for c in chunk.columns}
    data = {}
    n = len(chunk)
    missing = []
    for col in feature_cols:
        if col in cleaned:
            data[col] = pd.to_numeric(chunk[cleaned[col]], errors="coerce").to_numpy(dtype=float)
        elif col in fixed.index:
            try:
                val = float(fixed[col])
            except (TypeError, ValueError):
                val = float("nan")
            data[col] = np.full(n, val, dtype=float)
        else:
            missing.append(col)
    if missing:
        raise ScoringError(
            "Feature columns missing from the inputs: " + ", ".join(missing[:8])
            + (" ..." if len(missing) > 8 else ""), 8)
    X = pd.DataFrame(data, columns=feature_cols, index=chunk.index)

    # 3. rows with missing / non-finite descriptors. The tree models (XGBoost, LightGBM,
    #    CatBoost, scikit-learn RF >= 1.4) accept NaN, so by default such rows are scored as in
    #    SEARCH-ML v1; if a model refuses NaN, or nan_policy == "drop", they are skipped instead.
    values = X.to_numpy()
    finite = np.isfinite(values).all(axis=1)
    meta = base = None
    if not finite.all():
        if np.isinf(values).any():
            X = X.replace([np.inf, -np.inf], np.nan)
        if _STATE.get("nan_policy", "auto") != "drop":
            try:
                meta, base = ens.predict(X)
                counters["incomplete"] = int((~finite).sum())
            except ValueError:
                meta = base = None
        if meta is None:
            counters["invalid"] = int((~finite).sum())
            X = X.loc[finite]
            chunk = chunk.loc[finite]
    if X.empty:
        return _empty_result(id_col), counters

    # 4. predict
    if meta is None:
        meta, base = ens.predict(X)
    out = pd.DataFrame({
        "_pki": meta,
        id_col: chunk[id_col].astype(str).to_numpy(),
        META_COLS[0]: np.round(meta, 2),
        META_COLS[1]: np.round(KCAL_PER_PKI * meta, 2),
        META_COLS[2]: np.round(np.power(10.0, 9.0 - meta), 2),
    })
    for key in ("xgb", "cat", "rf", "lgbm"):
        out[MODEL_COLS[key]] = np.round(KCAL_PER_PKI * base[key], 2)
    counters["scored"] = int(len(out))
    return out, counters


def _empty_result(id_col: str) -> pd.DataFrame:
    return pd.DataFrame(columns=["_pki", id_col] + RESULT_COLS)


# --------------------------------------------------------------------------
# Results: external merge sort + summary statistics
# --------------------------------------------------------------------------

class SummaryStats:
    def __init__(self):
        self.hist = np.zeros(HIST_BINS, dtype=np.int64)
        self.n = 0
        self.sum = 0.0
        self.sumsq = 0.0
        self.min = math.inf
        self.max = -math.inf
        self.spread_sum = 0.0
        self.best: Optional[dict] = None

    def update(self, res: pd.DataFrame, id_col: str) -> None:
        if res.empty:
            return
        pki = res["_pki"].to_numpy(dtype=float)
        idx = np.clip(((pki - HIST_MIN) / HIST_STEP).astype(int), 0, HIST_BINS - 1)
        self.hist += np.bincount(idx, minlength=HIST_BINS)
        self.n += len(pki)
        self.sum += float(pki.sum())
        self.sumsq += float((pki * pki).sum())
        self.min = min(self.min, float(pki.min()))
        i = int(pki.argmax())
        if float(pki[i]) > self.max:
            self.max = float(pki[i])
            row = res.iloc[i]
            self.best = {"id": str(row[id_col]), "pki": round(float(pki[i]), 2),
                         "kcal": round(KCAL_PER_PKI * float(pki[i]), 2)}
        models = res[[MODEL_COLS[k] for k in MODEL_COLS]].to_numpy(dtype=float)
        self.spread_sum += float((models.max(axis=1) - models.min(axis=1)).sum())

    def as_dict(self) -> dict:
        if self.n == 0:
            return {"scored": 0, "histogram": None}
        mean = self.sum / self.n
        var = max(0.0, self.sumsq / self.n - mean * mean)
        return {
            "scored": self.n,
            "pki_mean": round(mean, 3),
            "pki_std": round(math.sqrt(var), 3),
            "pki_min": round(self.min, 3),
            "pki_max": round(self.max, 3),
            "kcal_best": round(KCAL_PER_PKI * self.max, 2),
            "mean_model_spread_kcal": round(self.spread_sum / self.n, 3),
            "best": self.best,
            "histogram": {
                "unit": "pKi",
                "start": HIST_MIN,
                "step": HIST_STEP,
                "counts": self.hist.tolist(),
            },
        }


class ResultCollector:
    """Accumulates scored rows and writes them globally sorted by predicted affinity."""

    def __init__(self, id_col: str, run_rows: int = 1_000_000, top: int = 0,
                 tmpdir: Optional[str] = None):
        self.id_col = id_col
        self.run_rows = max(500, int(run_rows))
        self.top = max(0, int(top or 0))
        self.buffer: List[pd.DataFrame] = []
        self.buffered = 0
        self.runs: List[str] = []
        self.tmpdir = tempfile.mkdtemp(prefix="searchml_runs_", dir=tmpdir)

    def add(self, res: pd.DataFrame) -> None:
        if res.empty:
            return
        self.buffer.append(res)
        self.buffered += len(res)
        if self.buffered >= self.run_rows:
            self._flush()

    def _sorted_buffer(self) -> pd.DataFrame:
        if not self.buffer:
            return _empty_result(self.id_col)
        df = pd.concat(self.buffer, ignore_index=True) if len(self.buffer) > 1 else self.buffer[0]
        df = df.sort_values("_pki", ascending=False, kind="mergesort")
        if self.top:
            df = df.head(self.top)
        return df

    def _flush(self) -> None:
        df = self._sorted_buffer()
        self.buffer, self.buffered = [], 0
        if df.empty:
            return
        path = os.path.join(self.tmpdir, f"run_{len(self.runs):05d}.csv")
        df.to_csv(path, index=False, header=False)
        self.runs.append(path)

    def finalize(self, output: str, preview: Optional[str] = None, preview_rows: int = 2000) -> int:
        header = [self.id_col] + RESULT_COLS
        written = 0
        try:
            if not self.runs:
                df = self._sorted_buffer().drop(columns=["_pki"])
                df.to_csv(output, index=False, columns=header)
                written = len(df)
                if preview:
                    df.head(preview_rows).to_csv(preview, index=False, columns=header)
                return written
            self._flush()
            handles = [open(p, "r", newline="", encoding="utf-8") for p in self.runs]
            try:
                readers = [csv.reader(h) for h in handles]
                merged = heapq.merge(*readers, key=lambda r: -float(r[0]))
                pv_fh = open(preview, "w", newline="", encoding="utf-8") if preview else None
                pv_w = csv.writer(pv_fh, lineterminator="\n") if pv_fh else None
                with open(output, "w", newline="", encoding="utf-8") as out_fh:
                    w = csv.writer(out_fh, lineterminator="\n")
                    w.writerow(header)
                    if pv_w:
                        pv_w.writerow(header)
                    for row in merged:
                        w.writerow(row[1:])
                        if pv_w and written < preview_rows:
                            pv_w.writerow(row[1:])
                        written += 1
                        if self.top and written >= self.top:
                            break
                if pv_fh:
                    pv_fh.close()
            finally:
                for h in handles:
                    h.close()
            return written
        finally:
            shutil.rmtree(self.tmpdir, ignore_errors=True)

    def cleanup(self) -> None:
        shutil.rmtree(self.tmpdir, ignore_errors=True)


# --------------------------------------------------------------------------
# Progress reporting
# --------------------------------------------------------------------------

class Progress:
    def __init__(self, path: Optional[str], total: int, estimated: bool, label: str,
                 percent_range: Tuple[float, float] = (0.0, 100.0), interval: float = 1.0):
        self.path = path
        self.total = int(total)
        self.estimated = estimated
        self.label = label
        self.lo, self.hi = percent_range
        self.interval = interval
        self.t0 = time.time()
        self.last = 0.0
        self.current = 0
        self.scored = 0

    def update(self, current: int, scored: int, force: bool = False, message: Optional[str] = None) -> None:
        self.current, self.scored = int(current), int(scored)
        if not self.path:
            return
        now = time.time()
        if not force and now - self.last < self.interval:
            return
        self.last = now
        elapsed = max(1e-6, now - self.t0)
        rate = self.current / elapsed
        frac = min(1.0, self.current / self.total) if self.total > 0 else 0.0
        eta = int((self.total - self.current) / rate) if rate > 0 and self.total > self.current else 0
        payload = {
            "stage": "scoring",
            "stage_label": self.label,
            "percent": round(self.lo + (self.hi - self.lo) * frac, 1),
            "stage_percent": round(100 * frac, 1),
            "current": self.current,
            "total": self.total,
            "total_is_estimate": self.estimated,
            "scored": self.scored,
            "rate_per_s": round(rate, 1),
            "eta_seconds": eta,
            "message": message or f"Scored {self.scored:,} of {self.current:,} rows read",
            "updated_at": int(now),
        }
        try:
            atomic_write_json(self.path, payload)
        except OSError:
            pass


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def run_screening(
    inputs: Sequence[str],
    fixed_csv: str,
    models_dir: str,
    output: str,
    role: str = "ligand",
    params_file: Optional[str] = None,
    order_file: Optional[str] = None,
    preview: Optional[str] = None,
    preview_rows: int = 2000,
    summary: Optional[str] = None,
    progress: Optional[str] = None,
    progress_range: Tuple[float, float] = (0.0, 100.0),
    workers: int = 1,
    threads: Optional[int] = None,
    top: int = 0,
    chunksize: int = 50_000,
    run_rows: int = 1_000_000,
    start_method: Optional[str] = None,
    total_rows: Optional[int] = None,
    nan_policy: str = "auto",
) -> dict:
    """Score ``inputs`` (molecules if role='ligand', proteins if role='protein') against the
    single row in ``fixed_csv``. Returns the summary dictionary."""
    t0 = time.time()
    paths = expand_inputs(inputs)

    if not os.path.isfile(fixed_csv):
        raise ScoringError(f"Feature file not found: {fixed_csv}", 2)
    fixed = drop_unnamed(pd.read_csv(fixed_csv))
    if fixed.empty:
        raise ScoringError(f"{fixed_csv} contains no rows.", 2)
    if len(fixed) != 1:
        log(f"Warning: {fixed_csv} has {len(fixed)} rows - using the first one.")
    fixed = fixed.iloc[[0]].reset_index(drop=True)

    # identifier column of the streamed side
    try:
        head = drop_unnamed(pd.read_csv(paths[0], nrows=5))
    except pd.errors.EmptyDataError:
        raise ScoringError(f"{paths[0]} is empty.", 2)
    if head.columns.empty:
        raise ScoringError(f"{paths[0]} has no columns.", 2)
    candidates = LIGAND_ID_CANDIDATES if role == "ligand" else PROTEIN_ID_CANDIDATES
    id_col = select_id_column(list(head.columns), candidates)

    filters = parse_parameter_file(params_file)
    if filters:
        log("Property filters: " + "; ".join(f"{k} {lo:g}-{hi:g}" for k, (lo, hi) in filters.items()))

    workers = max(1, int(workers or 1))
    if workers > 1:
        # keep each worker single-threaded; must be set before the model libraries load
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ.setdefault(var, str(max(1, int(threads or 1))))

    ensemble = Ensemble(models_dir, order_file)
    if ensemble.feature_columns is None:
        raise ScoringError("The scaler does not record feature names; cannot align features.", 8)

    if total_rows is not None and total_rows > 0:
        total, estimated = int(total_rows), False
    else:
        total, estimated = estimate_total_rows(paths)
    label = "Scoring molecules" if role == "ligand" else "Scoring protein targets"
    log(f"{label}: ~{total:,} rows in {len(paths)} file(s); workers={workers}, chunk={chunksize:,}")

    prog = Progress(progress, total, estimated, label, progress_range)
    prog.update(0, 0, force=True, message="Loading data")
    stats = SummaryStats()
    collector = ResultCollector(id_col, run_rows=run_rows, top=top,
                                tmpdir=os.path.dirname(os.path.abspath(output)) or None)
    counters = {"read": 0, "filtered": 0, "invalid": 0, "incomplete": 0, "scored": 0}

    _setup_state(ensemble, fixed, filters, id_col, nan_policy)

    def consume(result: Tuple[pd.DataFrame, Dict[str, int]]) -> None:
        res, c = result
        for k in counters:
            counters[k] += c[k]
        stats.update(res, id_col)
        collector.add(res)
        prog.update(counters["read"], counters["scored"])

    try:
        if workers == 1:
            ensemble.limit_threads(threads)
            for chunk in iter_chunks(paths, chunksize):
                consume(score_chunk(chunk))
        else:
            method = start_method or os.environ.get("SEARCHML_MP_START") or (
                "fork" if "fork" in mp.get_all_start_methods() else "spawn")
            ctx = mp.get_context(method)
            spawn_args = None
            if method != "fork":
                spawn_args = {"models_dir": models_dir, "order_file": order_file, "fixed": fixed,
                              "filters": filters, "id_col": id_col, "nan_policy": nan_policy}
            per_worker = max(1, int(threads or 1))
            pool = ctx.Pool(workers, initializer=_init_worker, initargs=(per_worker, spawn_args))
            try:
                inflight: deque = deque()
                max_inflight = workers * 2
                for chunk in iter_chunks(paths, chunksize):
                    inflight.append(pool.apply_async(score_chunk, (chunk,)))
                    while len(inflight) >= max_inflight:
                        consume(inflight.popleft().get())
                while inflight:
                    consume(inflight.popleft().get())
                pool.close()
                pool.join()
            except BaseException:
                pool.terminate()
                raise
    except BaseException:
        collector.cleanup()
        raise

    prog.update(counters["read"], counters["scored"], force=True, message="Sorting results")
    written = collector.finalize(output, preview, preview_rows)
    elapsed = time.time() - t0

    result = {
        "engine_version": __version__,
        "role": role,
        "id_column": id_col,
        "inputs": [os.path.basename(p) for p in paths],
        "rows_read": counters["read"],
        "filtered_out": counters["filtered"],
        "invalid_rows": counters["invalid"],
        "rows_with_missing_descriptors": counters["incomplete"],
        "rows_written": written,
        "top_limit": top or None,
        "filters": {k: list(v) for k, v in filters.items()},
        "model_order": ensemble.order,
        "elapsed_seconds": round(elapsed, 1),
        "workers": workers,
    }
    result.update(stats.as_dict())
    if summary:
        atomic_write_json(summary, result)
    prog.total = max(prog.total, counters["read"])
    prog.update(counters["read"], counters["scored"], force=True, message="Scoring finished")

    log(f"Rows read: {counters['read']:,} | filtered out: {counters['filtered']:,} | "
        f"invalid: {counters['invalid']:,} | scored: {counters['scored']:,} | written: {written:,}")
    if counters["incomplete"]:
        log(f"Note: {counters['incomplete']:,} rows had missing descriptors and were scored with the "
            "models' missing-value handling (use --drop-incomplete to skip them).")
    log(f"Results saved to '{output}' ({elapsed:.1f} s)")
    if counters["read"] > 0 and counters["scored"] == 0:
        if counters["filtered"] == counters["read"]:
            log("Note: every row was removed by the property filters - widen the ranges to get hits.")
        elif counters["invalid"]:
            log("Note: no row had a complete set of valid descriptors.")
    return result


def add_common_args(parser) -> None:
    """CLI options shared by screening.py and reverse_screening.py."""
    parser.add_argument("--models", default="./", help="Directory with the *.joblib models (default: ./)")
    parser.add_argument("--order", default=None, help=f"Path to {ORDER_FILE} (auto-detected if omitted)")
    parser.add_argument("--output", default="results.csv", help="Sorted output CSV (default: results.csv)")
    parser.add_argument("--preview", default=None, help="Also write the top rows to this CSV")
    parser.add_argument("--preview-rows", type=int, default=2000, help="Rows in the preview file (default 2000)")
    parser.add_argument("--summary", default=None, help="Write run statistics + histogram as JSON")
    parser.add_argument("--progress", default=None, help="Write live progress JSON to this path")
    parser.add_argument("--workers", type=int, default=int(os.environ.get("SEARCHML_WORKERS", "1")),
                        help="Parallel worker processes (default 1)")
    parser.add_argument("--threads", type=int, default=None,
                        help="Threads per model (default: all cores with 1 worker, 1 thread per worker otherwise)")
    parser.add_argument("--top", type=int, default=0, help="Keep only the N best rows (0 = keep all)")
    parser.add_argument("--chunksize", type=int, default=50_000, help="Rows per chunk (default 50000)")
    parser.add_argument("--total-rows", type=int, default=None, help="Known row count (for progress)")
    parser.add_argument("--drop-incomplete", action="store_true",
                        help="Skip rows with missing descriptors instead of scoring them")
    parser.add_argument("--version", action="version", version=f"SEARCH-ML scoring engine {__version__}")


def cli_main(run) -> int:
    try:
        run()
        return 0
    except ScoringError as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return exc.exit_code
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130
    except MemoryError:
        print("ERROR: out of memory while scoring. Use smaller --chunksize or fewer --workers.",
              file=sys.stderr)
        return 9
