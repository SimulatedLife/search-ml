import json
import os
import shutil
import sys
import tempfile
import unittest
import warnings

import pandas as pd

from helpers import REPO, SCRIPTS, make_fake_models, protein_row, reference_screen

sys.path.insert(0, SCRIPTS)
import scoring  # noqa: E402

FDA = os.path.join(REPO, "datasets", "fda", "dataset.csv")
DRUGBANK = os.path.join(REPO, "datasets", "drugbank", "dataset.csv")
PARAMS = os.path.join(REPO, "parameters", "parameter.txt")


class ScoringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("ignore")
        cls.tmp = tempfile.mkdtemp(prefix="searchml_test_")
        os.makedirs(os.path.join(cls.tmp, "models"))
        cls.models = make_fake_models(os.path.join(cls.tmp, "models"))
        cls.protein = protein_row(os.path.join(cls.tmp, "protein.csv"))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, name):
        return os.path.join(self.tmp, name)

    def run_screen(self, name, inputs=None, **kw):
        kw.setdefault("params_file", PARAMS)
        return scoring.run_screening(inputs or [FDA], self.protein, self.models, self.out(name), **kw)

    def test_matches_original_logic(self):
        summary = self.run_screen("a.csv")
        new = pd.read_csv(self.out("a.csv"))
        ref = reference_screen(self.models, FDA, self.protein, scoring.parse_parameter_file(PARAMS))
        self.assertEqual(len(new), len(ref))
        self.assertEqual(set(new["PDB_File"]), set(ref["PDB_File"]))
        merged = new.merge(ref, on="PDB_File", suffixes=("", "_ref"))
        for col in ("PBE (pKi/pKd/pIC50)", "PBE (kcal/mol)", "XGB (kcal/mol)", "RF (kcal/mol)"):
            self.assertEqual((merged[col] - merged[col + "_ref"]).abs().max(), 0.0, col)
        self.assertTrue(new["PBE (kcal/mol)"].is_monotonic_increasing)
        self.assertEqual(list(new.columns), ["PDB_File"] + scoring.RESULT_COLS)
        self.assertEqual(summary["rows_written"], len(new))
        self.assertEqual(summary["rows_read"], 2066)
        self.assertEqual(sum(summary["histogram"]["counts"]), len(new))

    def test_external_merge_parallel_and_top(self):
        self.run_screen("single.csv", [DRUGBANK], chunksize=1000)
        self.run_screen("multi.csv", [DRUGBANK], chunksize=700, run_rows=500, workers=2)
        a = open(self.out("single.csv")).read()
        b = open(self.out("multi.csv")).read()
        self.assertEqual(a, b, "external merge / parallel output differs from single pass")
        self.run_screen("top.csv", [DRUGBANK], chunksize=700, run_rows=500, top=250,
                        preview=self.out("prev.csv"), preview_rows=40)
        top = open(self.out("top.csv")).read().splitlines()
        self.assertEqual(top, a.splitlines()[:251])
        self.assertEqual(open(self.out("prev.csv")).read().splitlines(), a.splitlines()[:41])

    def test_no_filters_and_progress(self):
        prog = self.out("progress.json")
        s = self.run_screen("nf.csv", params_file=None, progress=prog, summary=self.out("s.json"))
        self.assertEqual(s["rows_written"], 2066)
        p = json.load(open(prog))
        self.assertEqual(p["current"], 2066)
        self.assertEqual(p["percent"], 100.0)
        self.assertTrue(json.load(open(self.out("s.json")))["best"]["id"])

    def test_filters_remove_everything(self):
        f = self.out("strict.txt")
        open(f, "w").write("MolWt: 1,2\n")
        s = self.run_screen("empty.csv", params_file=f)
        self.assertEqual(s["rows_written"], 0)
        self.assertEqual(open(self.out("empty.csv")).read().strip().split(",")[0], "PDB_File")

    def test_parameter_parsing(self):
        f = self.out("p.txt")
        open(f, "w").write("MolWt: 500,200\nTPSA: 10\n# comment\nbad line\nX: a,b\n")
        p = scoring.parse_parameter_file(f)
        self.assertEqual(p["MolWt"], (200.0, 500.0))
        self.assertEqual(p["TPSA"], (10.0, 10.0))
        self.assertNotIn("X", p)
        self.assertEqual(scoring.parse_parameter_file(self.out("missing.txt")), {})

    def test_model_order_default_matches_shipped_file(self):
        shipped = open(os.path.join(REPO, "parameters", "final_ensemble_model_order.txt")).read().split()
        self.assertEqual(scoring.DEFAULT_MODEL_ORDER, shipped)
        self.assertEqual(scoring.read_model_order(None), shipped)

    def test_missing_features_reported(self):
        bad = self.out("bad.csv")
        pd.read_csv(FDA, nrows=10).drop(columns=["TPSA"]).to_csv(bad, index=False)
        with self.assertRaises(scoring.ScoringError) as cm:
            self.run_screen("bad_out.csv", [bad], params_file=None)
        self.assertIn("TPSA", str(cm.exception))

    def test_invalid_rows_dropped(self):
        df = pd.read_csv(FDA, nrows=20)
        df["TPSA"] = df["TPSA"].astype(object)
        df.loc[3, "MolLogP"] = float("nan")
        df.loc[5, "TPSA"] = "oops"
        p = self.out("nan.csv")
        df.to_csv(p, index=False)
        # the stand-in linear models refuse NaN -> rows are skipped automatically
        s = self.run_screen("nan_out.csv", [p], params_file=None)
        self.assertEqual(s["invalid_rows"], 2)
        self.assertEqual(s["rows_written"], 18)
        s = self.run_screen("nan_out2.csv", [p], params_file=None, nan_policy="drop")
        self.assertEqual(s["rows_written"], 18)

    def test_reverse_role(self):
        lig = self.out("lig.csv")
        pd.read_csv(FDA, nrows=1).to_csv(lig, index=False)
        homo = self.out("homo.csv")
        pd.read_csv(os.path.join(REPO, "datasets", "HOMO", "protein.csv"), nrows=300).to_csv(homo, index=False)
        s = scoring.run_screening([homo], lig, self.models, self.out("rev.csv"), role="protein")
        out = pd.read_csv(self.out("rev.csv"))
        # HOMO contains entries without pocket descriptors; the stand-in models cannot score NaN
        self.assertEqual(len(out) + s["invalid_rows"], 300)
        self.assertTrue(out["PDB_File"].str.endswith(".pdb").all())
        self.assertEqual(s["role"], "protein")


if __name__ == "__main__":
    unittest.main()
