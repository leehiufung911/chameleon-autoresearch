import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys

sys.path.insert(0, "chameleon_local")

import numpy as np
from rdkit import Chem
import json

print("Starting iteration 19...", flush=True)

# Test file reading
protacs = []
with open("chameleon_local/user_protacs.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            protacs.append((parts[0], parts[1]))

print(f"Loaded {len(protacs)} PROTACs", flush=True)


# Simple spectral analysis
def analyze_mol(smiles, name):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()
    return {"name": name, "atoms": n_atoms, "bonds": n_bonds}


results = []
for name, smiles in protacs:
    res = analyze_mol(smiles, name)
    if res:
        results.append(res)
        print(f"{name}: {res['atoms']} atoms, {res['bonds']} bonds", flush=True)

print("Done with PROTACs", flush=True)

# Read benchmark (format: name smiles)
benchmark = []
with open("chameleon_local/benchmark.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            name = parts[0]
            smiles = parts[1]
            benchmark.append((name, smiles))

print(f"Loaded {len(benchmark)} benchmark molecules", flush=True)

print("SUCCESS: Iteration 19 test complete", flush=True)

# Write output file
with open("experiments/iter_19_output.txt", "w") as f:
    f.write(f"Loaded {len(protacs)} PROTACs\n")
    f.write(f"Loaded {len(benchmark)} benchmark molecules\n")
    f.write("SUCCESS: Complete\n")
    for r in results:
        f.write(f"{r['name']}: {r['atoms']} atoms, {r['bonds']} bonds\n")
