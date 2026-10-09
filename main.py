#!/usr/bin/env python3
"""Report residues within a distance cutoff of a target chain in a structure file.

Default query: lysine (K) residues within 4.0 A of chain C, the analysis used to
map protein-nucleic acid contacts in AlphaFold3-predicted complex models.
Distances are minimum atom-to-atom distances; hydrogen atoms are ignored.

Usage:
    python main.py --input your_model.cif [--chain C] [--resname K] [--cutoff 4.0]

Options:
    --input, -i PATH   input structure file (.cif/.mmcif or .pdb)
    --chain C          target chain id (default: C)
    --resname K        residue to search, one-letter or three-letter (default: LYS)
    --cutoff 4.0       distance cutoff in angstrom (default: 4.0)
    --top N            print only the N closest residues (0 = all)
    --models all       scan every model instead of the first one
    --csv PATH         also write the contacts to a CSV file

Requires Biopython:  pip install -r requirements.txt
"""
from __future__ import annotations

import argparse
import sys
import warnings

import numpy as np
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

CSV_HEADER = (
    "chain,residue,residue_number,insertion_code,min_distance_angstrom,"
    "query_atom,target_residue,target_atom"
)


def three_letter(name: str) -> str:
    """Normalize a one-letter or three-letter residue code to upper-case three letters.

    `protein_letters_1to3` maps to title case ("Lys"), while PDB/mmCIF residue names
    are upper case ("LYS"), so the mapped value needs `.upper()` as well.
    """
    token = name.strip().upper()
    if len(token) == 1:
        return protein_letters_1to3.get(token, token).upper()
    return token


def load_models(path: str, every_model: bool):
    """Parse `path` and return the models to scan (first model unless `every_model`)."""
    parser_cls = MMCIFParser if path.lower().endswith((".cif", ".mmcif")) else PDBParser
    structure = parser_cls(QUIET=True).get_structure("structure", path)
    models = list(structure.get_models())
    if not models:
        raise SystemExit(f"no model found in {path}")
    return models if every_model else models[:1]


def chain_summary(model) -> None:
    """Print atoms, residues and lysine counts for every chain of one model."""
    print("=== overview ===")
    for chain in model:
        residues = list(chain.get_residues())
        lysines = [r for r in residues if r.get_resname().strip() == "LYS"]
        ids = ",".join(str(r.id[1]) for r in lysines) or "none"
        print(f"chain {chain.id or '(blank)'}: {len(list(chain.get_atoms()))} atoms, "
              f"{len(residues)} residues, LYS={ids}")


def contacts(model, chain_id: str, resname: str, cutoff: float):
    """Return matching residues sorted by their minimum distance to `chain_id`.

    Each item is ((chain, resseq, icode), (distance, resname, query_atom, target_atom)).
    """
    target_atoms = [a for a in model.get_atoms()
                    if a.get_parent().get_parent().id == chain_id and a.element != "H"]
    if not target_atoms:
        raise SystemExit(f"chain {chain_id!r} has no atoms in this model")

    # KD-tree over the target chain; search() returns only atoms inside the cutoff,
    # so the closest returned atom is the residue's minimum distance to that chain.
    search = NeighborSearch(target_atoms)
    best: dict[tuple, tuple] = {}
    for residue in model.get_residues():
        chain = residue.get_parent().id
        if chain == chain_id or residue.get_resname().strip() != resname:
            continue
        for atom in residue:
            if atom.element == "H":
                continue
            for hit in search.search(atom.coord, cutoff):
                distance = float(np.linalg.norm(atom.coord - hit.coord))
                key = (chain, residue.id[1], residue.id[2].strip())
                if key not in best or distance < best[key][0]:
                    best[key] = (distance, residue.get_resname().strip(), atom, hit)
    return sorted(best.items(), key=lambda item: item[1][0])


def csv_rows(rows, chain_id: str):
    """Yield CSV lines for one model's contact rows."""
    for (chain, resseq, icode), (distance, resname, atom, hit) in rows:
        target = hit.get_parent()
        yield (f"{chain},{resname},{resseq},{icode or ''},{distance:.2f},"
               f"{atom.get_name()},{target.get_resname().strip()}{target.id[1]},{hit.get_name()}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Residues within a distance cutoff of a target chain (PDB/mmCIF).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("positional", nargs="?", help="input structure file (same as --input)")
    parser.add_argument("--input", "-i", dest="input_option", help="input structure file")
    parser.add_argument("--chain", default="C", help="target chain id (default: C)")
    parser.add_argument("--resname", default="LYS", help="residue name, 1- or 3-letter (default: LYS)")
    parser.add_argument("--cutoff", type=float, default=4.0, help="cutoff in angstrom (default: 4.0)")
    parser.add_argument("--top", type=int, default=0, help="print only the N closest (0 = all)")
    parser.add_argument("--models", choices=("first", "all"), default="first",
                        help="scan the first model or every model (default: first)")
    parser.add_argument("--csv", dest="csv_path", help="also write the contacts to this CSV file")
    args = parser.parse_args()

    path = args.input_option or args.positional
    if not path:
        parser.error("an input structure file is required (--input your_model.cif)")

    warnings.simplefilter("ignore", PDBConstructionWarning)
    residue_name = three_letter(args.resname)
    models = load_models(path, args.models == "all")

    hits = []
    csv_lines = [CSV_HEADER]
    for model in models:
        chain_summary(model)
        found = contacts(model, args.chain, residue_name, args.cutoff)
        prefix = f"model {model.id} " if args.models == "all" else ""
        print(f"\n{prefix}{residue_name} residues within <= {args.cutoff:g} A "
              f"of chain {args.chain}: {len(found)}")
        for (chain, resseq, icode), (distance, resname, atom, hit) in found:
            target = hit.get_parent()
            detail = (f"{atom.get_name()} vs {args.chain}:{target.get_resname().strip()}"
                      f"{target.id[1]}:{hit.get_name()}")
            print(f"{chain}:{resname}{resseq}{icode}  min={distance:.2f} A  ({detail})")
        hits.extend(found)
        csv_lines.extend(csv_rows(found, args.chain))

    if args.top > 0 and hits:
        print(f"\n=== {args.top} closest overall ===")
        for (chain, resseq, icode), (distance, resname, atom, hit) in sorted(hits, key=lambda i: i[1][0])[:args.top]:
            print(f"{chain}:{resname}{resseq}{icode}  {distance:.2f} A")

    if args.csv_path:
        with open(args.csv_path, "w", encoding="utf-8", newline="") as handle:
            handle.write("\n".join(csv_lines) + "\n")
        print(f"\nCSV written to {args.csv_path} ({len(hits)} contact(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
