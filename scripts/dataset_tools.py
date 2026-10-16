#!/usr/bin/env python
"""
Utilities for installing and maintaining SEARCH-ML libraries.

    # write count.txt (row count used for progress bars) for every installed library
    python dataset_tools.py count --catalog
    python dataset_tools.py count datasets/chembl

    # split a very large feature CSV into parts (faster restarts, parallel reading)
    python dataset_tools.py split /data/pubchem_features.csv datasets/pubchem --rows 2000000 [--gzip]

    # check that a library has every descriptor the models need
    python dataset_tools.py validate datasets/zinc --kind ligand
    python dataset_tools.py validate datasets/reverse/mus --kind protein

    # build a brand-new compound library from an SDF (3D, with hydrogens preferred)
    python dataset_tools.py featurize my_library.sdf datasets/mylib
    #   then add a line to config/databases.tsv:   mylib<TAB>MYLIB<TAB>My library<TAB>mylib<TAB>small
"""

import argparse
import csv
import gzip
import os
import sys

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)

LIGAND_FEATURES = [
    "LigVol", "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds",
    "NumAromaticRings", "FractionCSP3", "NumHeteroatoms", "HeavyAtomCount", "NumAtoms", "MolMR",
    "BertzCT", "NumCarbons", "Asphericity", "Eccentricity", "InertialShapeFactor", "NPR1", "NPR2",
    "PMI1", "PMI2", "PMI3", "RadiusOfGyration", "SpherocityIndex",
]
PROTEIN_FEATURES = [
    "protein_total_residues", "protein_molecular_weight", "protein_isoelectric_point",
    "protein_aromaticity", "protein_helix_fraction", "protein_turn_fraction", "protein_sheet_fraction",
    "pocket_num_residues", "pocket_heavy_atoms",
] + ["pocket_%s_count" % aa for aa in (
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO",
    "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR")] + [
    "pocket_polar_count", "pocket_nonpolar_count", "pocket_aromatic_count", "pocket_charged_pos_count",
    "pocket_charged_neg_count", "pocket_net_charge", "pocket_hydrophobicity_avg", "pocket_volume",
    "pocket_sasa",
]


def data_files(path):
    if os.path.isfile(path):
        return [path]
    for name in ("dataset.csv", "dataset.csv.gz", "protein.csv", "protein.csv.gz"):
        if os.path.isfile(os.path.join(path, name)):
            return [os.path.join(path, name)]
    parts = os.path.join(path, "parts") if os.path.isdir(os.path.join(path, "parts")) else path
    files = sorted(os.path.join(parts, f) for f in os.listdir(parts) if f.endswith((".csv", ".csv.gz")))
    if not files:
        raise SystemExit("No CSV data found in %s" % path)
    return files


def count_rows(path):
    total = 0
    for f in data_files(path):
        opener = gzip.open if f.endswith(".gz") else open
        with opener(f, "rb") as fh:
            n = 0
            for chunk in iter(lambda: fh.read(1 << 22), b""):
                n += chunk.count(b"\n")
        total += max(0, n - 1)
    return total


def cmd_count(a):
    targets = list(a.paths)
    if a.catalog:
        import pipeline
        for kind in ("databases", "libraries"):
            for e in pipeline.load_catalog(kind).values():
                p = e["path"] if os.path.isabs(e["path"]) else os.path.join(pipeline.datasets_dir(), e["path"])
                if os.path.exists(p):
                    targets.append(p)
    for p in targets:
        n = count_rows(p)
        dest = os.path.join(p if os.path.isdir(p) else os.path.dirname(p), "count.txt")
        try:
            with open(dest, "w") as fh:
                fh.write("%d\n" % n)
            print("%-60s %12s rows -> %s" % (p, format(n, ","), dest))
        except OSError as exc:
            print("%-60s %12s rows (could not write count.txt: %s)" % (p, format(n, ","), exc))
    return 0


def cmd_split(a):
    os.makedirs(os.path.join(a.outdir, "parts"), exist_ok=True)
    opener = gzip.open if a.input.endswith(".gz") else open
    part, rows, total = 0, 0, 0
    out = writer = None
    with opener(a.input, "rt", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        for row in reader:
            if out is None or rows >= a.rows:
                if out:
                    out.close()
                part += 1
                name = os.path.join(a.outdir, "parts", "part_%04d.csv%s" % (part, ".gz" if a.gzip else ""))
                out = gzip.open(name, "wt", newline="") if a.gzip else open(name, "w", newline="")
                writer = csv.writer(out, lineterminator="\n")
                writer.writerow(header)
                rows = 0
            writer.writerow(row)
            rows += 1
            total += 1
    if out:
        out.close()
    with open(os.path.join(a.outdir, "count.txt"), "w") as fh:
        fh.write("%d\n" % total)
    print("Wrote %d parts (%s rows) to %s/parts" % (part, format(total, ","), a.outdir))
    return 0


def cmd_validate(a):
    need = LIGAND_FEATURES if a.kind == "ligand" else PROTEIN_FEATURES
    ok = True
    for f in data_files(a.path)[:3]:
        opener = gzip.open if f.endswith(".gz") else open
        with opener(f, "rt", newline="") as fh:
            header = next(csv.reader(fh))
        missing = [c for c in need if c not in header]
        if missing:
            ok = False
            print("%s: MISSING %d column(s): %s" % (f, len(missing), ", ".join(missing)))
        else:
            print("%s: OK (identifier column '%s', %d columns)" % (f, header[0], len(header)))
    return 0 if ok else 1


def cmd_featurize(a):
    from ligand_features import calculate_features_sequentially
    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "dataset.csv")
    processed, failed = calculate_features_sequentially(a.sdf, out, prepare=not a.no_prepare)
    with open(os.path.join(a.outdir, "count.txt"), "w") as fh:
        fh.write("%d\n" % processed)
    print("Library written to %s (%d molecules, %d skipped)" % (out, processed, failed))
    return 0 if processed else 6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("count", help="write count.txt")
    c.add_argument("paths", nargs="*")
    c.add_argument("--catalog", action="store_true", help="all installed libraries from config/*.tsv")
    c.set_defaults(func=cmd_count)
    s = sub.add_parser("split", help="split a large CSV into parts/")
    s.add_argument("input")
    s.add_argument("outdir")
    s.add_argument("--rows", type=int, default=2_000_000)
    s.add_argument("--gzip", action="store_true")
    s.set_defaults(func=cmd_split)
    v = sub.add_parser("validate", help="check descriptor columns")
    v.add_argument("path")
    v.add_argument("--kind", choices=["ligand", "protein"], default="ligand")
    v.set_defaults(func=cmd_validate)
    f = sub.add_parser("featurize", help="build a library from an SDF")
    f.add_argument("sdf")
    f.add_argument("outdir")
    f.add_argument("--no-prepare", action="store_true")
    f.set_defaults(func=cmd_featurize)
    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
