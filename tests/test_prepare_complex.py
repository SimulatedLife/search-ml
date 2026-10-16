import os
import shutil
import sys
import tempfile
import unittest

from helpers import REPO, SCRIPTS

sys.path.insert(0, SCRIPTS)
import prepare_complex as pc  # noqa: E402

EXAMPLE = os.path.join(REPO, "example", "9kte.pdb")


class PrepareComplexTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_selector_formats(self):
        self.assertEqual(pc.parse_selector("9UM"), ("9UM", None, None, None))
        self.assertEqual(pc.parse_selector("9um:A:413"), ("9UM", "A", "413", ""))
        self.assertEqual(pc.parse_selector("VWW:_:210"), ("VWW", "", "210", ""))
        self.assertEqual(pc.parse_selector("VWW A 210"), ("VWW", "A", "210", None))
        self.assertEqual(pc.parse_selector("VWW A 210")[:3], ("VWW", "A", "210"))
        self.assertEqual(pc.parse_selector("BLA A1301")[:3], ("BLA", "A", "1301"))
        with self.assertRaises(pc.PrepareError):
            pc.parse_selector("")
        with self.assertRaises(pc.PrepareError):
            pc.parse_selector("bad;rm")

    def test_list_ligands(self):
        groups = pc.list_ligands(pc.first_altloc_only(pc.read_pdb_lines(EXAMPLE)))
        selectors = [g["selector"] for g in groups]
        self.assertIn("9UM:A:413", selectors)
        self.assertNotIn("HOH:A:501", selectors)
        ions = [g for g in groups if g["kind"] == "ion"]
        self.assertTrue(any(g["resname"] == "MG" for g in ions))

    def test_prepare_writes_protein_and_one_ligand(self):
        out = os.path.join(self.tmp, "PROT.pdb")
        info = pc.prepare(EXAMPLE, "9UM", out)
        self.assertEqual(info["ligand"]["resname"], "9UM")
        self.assertGreaterEqual(info["protein_residues"], 25)
        text = open(out).read().splitlines()
        lig = [ln for ln in text if ln[17:20] == "9UM"]
        self.assertEqual(len(lig), info["ligand"]["atoms"])
        self.assertFalse(any(ln[17:20] == "HOH" for ln in text))
        self.assertTrue(all(ln[16] == " " for ln in text if ln.startswith(("ATOM", "HETATM"))))

    def test_ligand_not_found_lists_alternatives(self):
        with self.assertRaises(pc.PrepareError) as cm:
            pc.prepare(EXAMPLE, "XYZ", os.path.join(self.tmp, "x.pdb"))
        self.assertEqual(cm.exception.exit_code, 13)
        self.assertIn("9UM:A:413", str(cm.exception))

    def test_short_protein_rejected(self):
        lines = [ln for ln in open(EXAMPLE) if ln.startswith(("ATOM", "HETATM"))]
        short, seen = [], set()
        for ln in lines:
            if ln.startswith("ATOM"):
                seen.add(ln[22:26])
                if len(seen) > 10:
                    continue
            short.append(ln)
        p = os.path.join(self.tmp, "short.pdb")
        open(p, "w").writelines(short)
        with self.assertRaises(pc.PrepareError) as cm:
            pc.prepare(p, "9UM", os.path.join(self.tmp, "o.pdb"))
        self.assertEqual(cm.exception.exit_code, 3)

    def test_altloc_first_only(self):
        a = "ATOM      1  CA ALYS A   1      11.000  11.000  11.000  0.50 10.00           C\n"
        b = "ATOM      2  CA BLYS A   1      12.000  12.000  12.000  0.50 10.00           C\n"
        kept = pc.first_altloc_only([a.rstrip("\n").ljust(80), b.rstrip("\n").ljust(80)])
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0][16], " ")


if __name__ == "__main__":
    unittest.main()
