"""Tests for hpc/bin/searchml-hpc with the 'local' scheduler in simulation mode."""
import json
import os
import shutil
import subprocess
import tempfile
import time
import unittest

from helpers import REPO

HPC = os.path.join(REPO, "hpc", "bin", "searchml-hpc")
GATEWAY = os.path.join(REPO, "hpc", "bin", "searchml-gateway")
EXAMPLE = os.path.join(REPO, "example", "9kte.pdb")


class DispatcherTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.jobs = os.path.join(cls.tmp, "jobs")
        cls.conf = os.path.join(cls.tmp, "hpc.conf")
        with open(cls.conf, "w") as fh:
            fh.write("SEARCH_ML_HOME=%s\nJOBS_DIR=%s\nPYTHON=%s\nSCHEDULER=local\nMAX_LOCAL_JOBS=2\n"
                     % (REPO, cls.jobs, shutil.which("python3")))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def call(self, *args, stdin=None, env=None, exe=HPC):
        e = os.environ.copy()
        e["SEARCHML_HPC_CONF"] = self.conf
        e.update(env or {})
        p = subprocess.run([exe] + list(args), input=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=e)
        return p.returncode, p.stdout

    def j(self, *args, **kw):
        rc, out = self.call(*args, **kw)
        return rc, json.loads(out.decode() or "{}")

    def submit(self, job_id, spec):
        self.assertEqual(self.j("prepare", job_id)[0], 0)
        self.assertEqual(self.j("put", job_id, "job.json", stdin=json.dumps(spec).encode())[0], 0)
        self.assertEqual(self.j("put", job_id, "receptor.pdb", stdin=open(EXAMPLE, "rb").read())[0], 0)
        rc, res = self.j("submit", job_id)
        self.assertEqual(rc, 0, res)
        return res

    def wait(self, job_id, timeout=40):
        t0 = time.time()
        while time.time() - t0 < timeout:
            rc, st = self.j("status", job_id)
            if st.get("state") in ("COMPLETED", "FAILED", "CANCELLED"):
                return st
            time.sleep(0.5)
        self.fail("job did not finish")

    def spec(self, **opts):
        o = {"simulate": True, "simulate_seconds": 1.5}
        o.update(opts)
        return {"version": 2, "mode": "forward",
                "forward": {"receptor": "receptor.pdb", "ligand": {"resname": "9UM", "chain": "A", "resseq": "413"},
                            "source": "database", "database": "fda", "filters": {}}, "options": o}

    def test_ping(self):
        rc, d = self.j("ping")
        self.assertEqual(rc, 0)
        self.assertTrue(d["ok"])
        self.assertTrue(d["databases"]["fda"]["available"])
        self.assertFalse(d["databases"]["pubchem"]["available"])

    def test_full_job_cycle(self):
        self.submit("test-job-0001", self.spec())
        st = self.wait("test-job-0001")
        self.assertEqual(st["state"], "COMPLETED", st)
        self.assertGreater(st["results"]["preview"], 100)
        rc, data = self.call("get", "test-job-0001", "results_preview.csv")
        self.assertEqual(rc, 0)
        self.assertTrue(data.startswith(b"PDB_File,"))

    def test_failure_reported(self):
        self.submit("test-job-0002", self.spec(simulate_fail="features"))
        st = self.wait("test-job-0002")
        self.assertEqual(st["state"], "FAILED")
        self.assertIn("Simulated failure", st["error"]["message"])

    def test_cancel(self):
        self.submit("test-job-0003", self.spec(simulate_seconds=20))
        time.sleep(1.5)
        rc, res = self.j("cancel", "test-job-0003")
        self.assertEqual(rc, 0)
        self.assertEqual(self.wait("test-job-0003")["state"], "CANCELLED")

    def test_input_validation(self):
        self.assertEqual(self.j("prepare", "../etc")[0], 1)
        self.j("prepare", "test-job-0004")
        rc, res = self.j("put", "test-job-0004", "shell.php", stdin=b"<?php")
        self.assertEqual(rc, 1)
        rc, res = self.j("get", "test-job-0004", "../../hpc.conf")
        self.assertEqual(rc, 1)
        rc, res = self.j("submit", "test-job-0004")
        self.assertEqual(rc, 1)
        self.assertIn("job.json", res["error"])

    def test_gateway(self):
        for cmd, ok in [("searchml-hpc ping", True), ("bash", False), ("", False),
                        ("searchml-hpc selftest", False), ("scp -f /etc/passwd", False),
                        ("searchml-hpc ping; id", False)]:
            rc, out = self.call(env={"SSH_ORIGINAL_COMMAND": cmd}, exe=GATEWAY)
            self.assertEqual(rc == 0, ok, (cmd, out))


if __name__ == "__main__":
    unittest.main()
