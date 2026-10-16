"""Shared test helpers: tiny scikit-learn stand-ins for the four base models.

The real scaler and meta-model from models/ are used so feature alignment and
the meta-model stacking order are exercised exactly as in production; XGBoost,
CatBoost, LightGBM and the 830 MB Random Forest are replaced by small sklearn
regressors trained on random data.
"""
import os
import shutil
import tempfile
import warnings

import joblib
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
SCRIPTS = os.path.join(REPO, "scripts")


def make_fake_models(dest=None):
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import LinearRegression, Ridge
    from sklearn.tree import DecisionTreeRegressor

    dest = dest or tempfile.mkdtemp(prefix="searchml_models_")
    src = os.path.join(REPO, "models")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scaler = joblib.load(os.path.join(src, "scaler_final.joblib"))
    cols = list(scaler.feature_names_in_)
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.normal(size=(400, len(cols))), columns=cols)
    y = 5 + 0.3 * X.iloc[:, :5].sum(axis=1)
    joblib.dump(LinearRegression().fit(X, y), os.path.join(dest, "final_xgb_model.joblib"))
    joblib.dump(DecisionTreeRegressor(max_depth=6, random_state=0).fit(X, y),
                os.path.join(dest, "final_cat_model.joblib"))
    joblib.dump(Ridge(alpha=3).fit(X, y), os.path.join(dest, "final_lgbm_model.joblib"))
    joblib.dump(RandomForestRegressor(n_estimators=10, random_state=0).fit(X, y),
                os.path.join(dest, "final_rf_model.joblib"))
    for f in ("scaler_final.joblib", "final_meta_model.joblib"):
        shutil.copy(os.path.join(src, f), dest)
    shutil.copy(os.path.join(REPO, "parameters", "final_ensemble_model_order.txt"), dest)
    return dest


def protein_row(dest_csv, index=1):
    df = pd.read_csv(os.path.join(REPO, "datasets", "HOMO", "protein.csv"), nrows=index + 1)
    df.iloc[[index]].to_csv(dest_csv, index=False)
    return dest_csv


def reference_screen(models_dir, molecules_csv, protein_csv, params):
    """Re-implementation of the original (v1) screening.py logic, used as the oracle."""
    import re
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scaler = joblib.load(os.path.join(models_dir, "scaler_final.joblib"))
        m = {k: joblib.load(os.path.join(models_dir, "final_%s_model.joblib" % k))
             for k in ("xgb", "cat", "rf", "lgbm")}
        meta = joblib.load(os.path.join(models_dir, "final_meta_model.joblib"))
    order = open(os.path.join(models_dir, "final_ensemble_model_order.txt")).read().split()
    mol = pd.read_csv(molecules_csv)
    prot = pd.read_csv(protein_csv).iloc[[0]]
    comb = pd.concat([mol.reset_index(drop=True),
                      pd.concat([prot] * len(mol), ignore_index=True)], axis=1)
    comb.columns = [re.sub(r"[\{\}\[\]:\"]", "", re.sub(r"[,\(\)\-\s]", "_", c)) for c in comb.columns]
    cols = list(scaler.feature_names_in_)
    X = comb.loc[:, ~comb.columns.duplicated()][cols]
    Xs = pd.DataFrame(scaler.transform(X), columns=cols)
    base = {k: v.predict(Xs) for k, v in m.items()}
    pki = meta.predict(np.column_stack([base[k] for k in order]))
    res = mol.copy()
    res["PBE (pKi/pKd/pIC50)"] = np.round(pki, 2)
    res["PBE (kcal/mol)"] = np.round(-1.364 * pki, 2)
    for k, name in (("xgb", "XGB"), ("cat", "CAT"), ("rf", "RF"), ("lgbm", "LGBM")):
        res["%s (kcal/mol)" % name] = np.round(-1.364 * base[k], 2)
    for p, (lo, hi) in params.items():
        if p in res.columns:
            res = res[(res[p] >= lo) & (res[p] <= hi)]
    return res
