"""End-to-end orchestration tests in simulation mode (no RDKit / AmberTools / models needed)."""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

from helpers import REPO, SCRIPTS

PIPELINE = os.path.join(SCRIPTS, "pipeline.py")
EXAMPLE = os.path.join(REPO, "example", "9kte.pdb")


def run(args, cwd, env_extra=None):
    env = os.environ.copy()
    env["SEARCH_ML_HOME"] = REPO
    env.update(env_extra or {})
    p = subprocess.run([sys.executable, PIPELINE] + args, cwd=cwd, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    return p.returncode, p.stdout


class PipelineSimulationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        shutil.copy(EXAMPLE, self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def job(self, name):
        return os.path.join(self.tmp, name)

    def test_forward_database(self):
        rc, out = run(["forward", "--pdb", "9kte.pdb", "--ligand", "9UM:A:413", "--mode", "FDA",
                       "--job-dir", "j1"], self.tmp, {"SEARCHML_SIMULATE": "1"})
        self.assertEqual(rc, 0, out)
        st = json.load(open(os.path.join(self.job("j1"), "state.json")))
        self.assertEqual(st["state"], "COMPLETED")
        self.assertTrue(os.path.exists(os.path.join(self.job("j1"), "COMPLETED")))
        self.assertFalse(os.path.exists(os.path.join(self.job("j1"), "RUNNING")))
        summ = json.load(open(os.path.join(self.job("j1"), "summary.json")))
        self.assertEqual(summ["rows_written"], 2066)
        head = open(os.path.join(self.job("j1"), "results.csv")).readline()
        self.assertTrue(head.startswith("PDB_File,PBE (pKi/pKd/pIC50)"))

    def test_reverse_smiles(self):
        rc, out = run(["reverse", "--smiles", "CC(=O)Oc1ccccc1C(=O)O", "--job-dir", "j2", "--top", "50"],
                      self.tmp, {"SEARCHML_SIMULATE": "1"})
        self.assertEqual(rc, 0, out)
        rows = open(os.path.join(self.job("j2"), "results.csv")).read().splitlines()
        self.assertEqual(len(rows), 51)

    def test_cli_exit_codes(self):
        self.assertEqual(run(["forward", "--pdb", "nope.pdb", "--ligand", "X", "--mode", "DB"], self.tmp)[0], 2)
        self.assertEqual(run(["forward", "--pdb", "9kte.pdb", "--ligand", "9UM", "--mode", "BOGUS"], self.tmp)[0], 6)
        self.assertEqual(run(["forward", "--pdb", "9kte.pdb", "--ligand", "9UM", "--mode", "SINGLE"], self.tmp)[0], 4)
        self.assertEqual(run(["forward", "--pdb", "9kte.pdb", "--ligand", "9UM", "--mode", "CUSTOM"], self.tmp)[0], 5)

    def test_simulated_failure_writes_error_log(self):
        d = self.job("jf")
        os.makedirs(d)
        shutil.copy(EXAMPLE, os.path.join(d, "receptor.pdb"))
        json.dump({"version": 2, "mode": "forward",
                   "forward": {"receptor": "receptor.pdb", "ligand": "9UM", "source": "database", "database": "fda"},
                   "options": {"simulate": True, "simulate_seconds": 1, "simulate_fail": "tleap"}},
                  open(os.path.join(d, "job.json"), "w"))
        rc, out = run(["run", d, "--quiet"], self.tmp)
        self.assertNotEqual(rc, 0)
        st = json.load(open(os.path.join(d, "state.json")))
        self.assertEqual(st["state"], "FAILED")
        self.assertEqual(st["stage"], "tleap")
        self.assertIn("Simulated failure", open(os.path.join(d, "error.log")).read())

    def test_bad_spec_rejected(self):
        d = self.job("jb")
        os.makedirs(d)
        json.dump({"version": 2, "mode": "sideways"}, open(os.path.join(d, "job.json"), "w"))
        rc, _ = run(["run", d, "--quiet", "--simulate"], self.tmp)
        self.assertEqual(rc, 6)

    def test_catalog(self):
        rc, out = run(["catalog", "--json"], self.tmp)
        self.assertEqual(rc, 0)
        cat = json.loads(out)
        keys = {d["key"]: d for d in cat["databases"]}
        self.assertTrue(keys["fda"]["available"])
        self.assertTrue({d["key"]: d for d in cat["libraries"]}["homo"]["available"])


if __name__ == "__main__":
    unittest.main()
