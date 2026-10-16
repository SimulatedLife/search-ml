#!/usr/bin/env python
"""
Prepare a protein-ligand complex for SEARCH-ML.

* Keeps the standard / AMBER amino-acid residues (ATOM records) and the one
  reference ligand that defines the binding pocket.
* Keeps only the first alternate location of every atom (like pdb4amber).
* Validates that the protein has at least 25 residues.
* Writes PROT.pdb, the input for AmberTools tleap (scripts/leap.cmd).

Usage
-----
    python prepare_complex.py --list complex.pdb            # show candidate ligands
    python prepare_complex.py complex.pdb VWW               # first VWW residue
    python prepare_complex.py complex.pdb VWW:A:210         # exact residue (RES:CHAIN:RESNUM)
    python prepare_complex.py complex.pdb "VWW A 210"       # PDB columns 18-26 form
    python prepare_complex.py complex.pdb VWW -o PROT.pdb --json prepare.json

Exit codes: 0 ok, 2 file not found / unreadable, 3 protein shorter than 25
residues, 13 reference ligand not found.
"""

import argparse
import json
import os
import re
import sys

# Residue names accepted as protein (standard + AMBER protonation variants)
PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY",
    "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
}
WATER = {"HOH", "WAT", "SOL", "H2O", "DOD", "TIP", "TIP3", "SPC"}
MIN_PROTEIN_RESIDUES = 25


class PrepareError(Exception):
    def __init__(self, message, exit_code):
        super().__init__(message)
        self.exit_code = exit_code


def parse_selector(text):
    """Parse a ligand selector into (resname, chain or None, resseq or None, icode or None)."""
    if text is None or not str(text).strip():
        raise PrepareError("No reference ligand given.", 13)
    raw = str(text)
    s = raw.strip()
    if ":" in s:
        parts = [p.strip() for p in s.split(":")]
        resname = parts[0].upper()
        chain = parts[1] if len(parts) > 1 else None
        if chain in ("", "_", "-", "N/A"):
            chain = "" if len(parts) > 1 else None
        resseq, icode = None, None
        if len(parts) > 2 and parts[2]:
            m = re.fullmatch(r"(-?\d+)([A-Za-z]?)", parts[2])
            if not m:
                raise PrepareError(f"Invalid residue number in ligand selector '{raw}'.", 13)
            resseq, icode = m.group(1), m.group(2).upper()
    elif len(raw) == 9 and raw[3] == " ":
        # exact PDB columns 18-26: RES(3) blank(1) CHAIN(1) RESSEQ(4)
        resname = raw[0:3].strip().upper()
        chain = raw[4].strip()
        resseq = raw[5:9].strip() or None
        icode = None
    else:
        tokens = s.split()
        resname = tokens[0].upper()
        chain, resseq, icode = None, None, None
        if len(tokens) == 2:
            m = re.fullmatch(r"([A-Za-z0-9]?)(-?\d+)([A-Za-z]?)", tokens[1])
            if m:
                chain, resseq, icode = m.group(1), m.group(2), m.group(3).upper()
        elif len(tokens) >= 3:
            chain = tokens[1]
            m = re.fullmatch(r"(-?\d+)([A-Za-z]?)", tokens[2])
            if m:
                resseq, icode = m.group(1), m.group(2).upper()
    if not re.fullmatch(r"[A-Z0-9+\-]{1,4}", resname):
        raise PrepareError(f"Invalid ligand residue name '{resname}'.", 13)
    return resname, chain, resseq, icode


def residue_key(line):
    return (line[17:20].strip().upper(), line[21:22].strip(), line[22:26].strip(), line[26:27].strip().upper())


def read_pdb_lines(path):
    if not os.path.isfile(path):
        raise PrepareError(f"PDB file '{path}' not found.", 2)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\r\n") for ln in fh]
    atoms = [ln.ljust(80) for ln in lines if ln.startswith(("ATOM  ", "HETATM", "ATOM ")) and len(ln) >= 54]
    if not atoms:
        raise PrepareError(f"'{path}' contains no ATOM/HETATM records.", 2)
    # keep the first model only
    if any(ln.startswith("MODEL") for ln in lines):
        first, in_model, seen_model = [], False, False
        for ln in lines:
            if ln.startswith("MODEL"):
                if seen_model:
                    break
                in_model, seen_model = True, True
            elif ln.startswith("ENDMDL"):
                break
            elif in_model and ln.startswith(("ATOM", "HETATM")) and len(ln) >= 54:
                first.append(ln.ljust(80))
        if first:
            atoms = first
    return atoms


def first_altloc_only(lines):
    """Keep atoms with a blank altloc or the first altloc seen for that residue."""
    chosen = {}
    out = []
    for ln in lines:
        alt = ln[16]
        key = residue_key(ln)
        if alt != " ":
            first = chosen.setdefault(key, alt)
            if alt != first:
                continue
            ln = ln[:16] + " " + ln[17:]
        out.append(ln)
    return out


def list_ligands(lines):
    groups = {}
    order = []
    for ln in lines:
        resname = ln[17:20].strip().upper()
        if not resname or resname in WATER:
            continue
        if ln.startswith("ATOM") and resname in PROTEIN_RESIDUES:
            continue
        key = residue_key(ln)
        if key not in groups:
            groups[key] = {"resname": key[0], "chain": key[1], "resseq": key[2], "icode": key[3],
                           "atoms": 0, "record": ln[:6].strip()}
            order.append(key)
        groups[key]["atoms"] += 1
    result = []
    for k in order:
        g = groups[k]
        g["selector"] = f"{g['resname']}:{g['chain'] or '_'}:{g['resseq']}{g['icode']}"
        g["kind"] = "ion" if g["atoms"] == 1 else "ligand"
        result.append(g)
    return result


def prepare(pdb_path, selector, out_path="PROT.pdb"):
    lines = first_altloc_only(read_pdb_lines(pdb_path))
    resname, chain, resseq, icode = parse_selector(selector)

    protein = [ln for ln in lines if ln.startswith("ATOM") and ln[17:20].strip().upper() in PROTEIN_RESIDUES]
    ca = {(ln[17:20], ln[21], ln[22:27]) for ln in protein if ln[12:16].strip() == "CA"}
    n_res = len(ca)
    if n_res < MIN_PROTEIN_RESIDUES:
        raise PrepareError(
            f"Protein size must be >= {MIN_PROTEIN_RESIDUES} amino acids. Found: {n_res}.", 3)

    def matches(ln):
        k = residue_key(ln)
        if k[0] != resname:
            return False
        if chain is not None and k[1] != chain:
            return False
        if resseq is not None and k[2] != resseq:
            return False
        if icode and k[3] != icode:
            return False
        return True

    candidates = [ln for ln in lines if ln.startswith(("ATOM", "HETATM")) and matches(ln)
                  and not (ln.startswith("ATOM") and resname in PROTEIN_RESIDUES and chain is None)]
    if not candidates:
        available = ", ".join(g["selector"] for g in list_ligands(lines)[:15]) or "none"
        raise PrepareError(
            f"Reference ligand '{selector}' not found in {os.path.basename(pdb_path)}. "
            f"Het groups present: {available}.", 13)
    instances = []
    for ln in candidates:
        k = residue_key(ln)
        if k not in instances:
            instances.append(k)
    selected = instances[0]
    warnings = []
    if len(instances) > 1:
        others = ", ".join(f"{k[0]}:{k[1] or '_'}:{k[2]}{k[3]}" for k in instances[1:6])
        warnings.append(
            f"{len(instances)} copies of {resname} found; using {selected[0]}:{selected[1] or '_'}:"
            f"{selected[2]}{selected[3]} (others: {others}). Use RES:CHAIN:RESNUM to choose another.")
    ligand = [ln for ln in candidates if residue_key(ln) == selected]

    # write protein (TER between chains) + ligand
    out = []
    prev_chain = None
    for ln in protein:
        if prev_chain is not None and ln[21] != prev_chain:
            out.append("TER")
        prev_chain = ln[21]
        out.append(ln.rstrip())
    out.append("TER")
    out.extend(ln.rstrip() for ln in ligand)
    out.append("END")
    with open(out_path, "w") as fh:
        fh.write("\n".join(out) + "\n")

    chains = sorted({ln[21].strip() or "_" for ln in protein})
    return {
        "input": os.path.basename(pdb_path),
        "output": out_path,
        "protein_residues": n_res,
        "protein_atoms": len(protein),
        "chains": chains,
        "ligand": {"resname": selected[0], "chain": selected[1], "resseq": selected[2],
                   "icode": selected[3], "atoms": len(ligand)},
        "warnings": warnings,
    }


def main():
    ap = argparse.ArgumentParser(description="Prepare a protein-ligand complex for SEARCH-ML (writes PROT.pdb).")
    ap.add_argument("pdb", help="Input PDB file")
    ap.add_argument("ligand", nargs="?", help="Reference ligand: RES, RES:CHAIN:RESNUM or 'RES C NNNN'")
    ap.add_argument("-o", "--output", default="PROT.pdb", help="Output PDB (default PROT.pdb)")
    ap.add_argument("--list", action="store_true", help="List candidate ligands / het groups and exit")
    ap.add_argument("--json", help="Write a JSON summary to this file")
    a = ap.parse_args()
    try:
        if a.list:
            groups = list_ligands(first_altloc_only(read_pdb_lines(a.pdb)))
            if not groups:
                print("No het groups (ligands/ions) found.")
            for g in groups:
                print(f"{g['selector']:<16} {g['kind']:<7} {g['atoms']:>4} atoms  ({g['record']})")
            if a.json:
                with open(a.json, "w") as fh:
                    json.dump(groups, fh, indent=1)
            return 0
        if not a.ligand:
            ap.error("the ligand argument is required (or use --list)")
        info = prepare(a.pdb, a.ligand, a.output)
        lig = info["ligand"]
        print(f"Protein length ({info['protein_residues']}) is OK. Chains: {', '.join(info['chains'])}")
        print(f"Reference ligand: {lig['resname']} chain {lig['chain'] or '-'} residue "
              f"{lig['resseq']}{lig['icode']} ({lig['atoms']} atoms)")
        for w in info["warnings"]:
            print("Warning: " + w)
        if a.json:
            with open(a.json, "w") as fh:
                json.dump(info, fh, indent=1)
        print(f"Wrote {a.output}")
        return 0
    except PrepareError as exc:
        print(f"Error: {exc}")
        return exc.exit_code


if __name__ == "__main__":
    sys.exit(main())
