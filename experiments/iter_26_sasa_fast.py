"""
Iteration 26: SASA Fluctuation Analysis (Fast Version)
Test SASA variance across small conformer ensembles for PROTACs
"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

print("=" * 70, file=sys.stderr)
print("ITERATION 26: SASA FLUCTUATION ENTROPY", file=sys.stderr)
print("=" * 70, file=sys.stderr)

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
import json

# Load PROTACs
df = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t")

results = []
output_lines = []

for _, row in df.iterrows():
    name = row["Name"]
    label = row["Label"]
    smiles = row["SMILES"]

    line = f"\n{name} ({label})"
    print(line, file=sys.stderr)
    output_lines.append(line)

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue

        mol = Chem.AddHs(mol)

        # Generate small ensemble
        AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=42)

        # Minimize
        for i in range(mol.GetNumConformers()):
            AllChem.MMFFOptimizeMolecule(mol, confId=i)

        # Compute SASA for each conformer
        radii = rdFreeSASA.classifyAtoms(mol)
        sasa_vals = []
        for i in range(mol.GetNumConformers()):
            sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=i)
            sasa_vals.append(sasa)

        sasa_vals = np.array(sasa_vals)

        stats = {
            "name": name,
            "label": label,
            "sasa_mean": float(np.mean(sasa_vals)),
            "sasa_std": float(np.std(sasa_vals)),
            "sasa_range": float(np.max(sasa_vals) - np.min(sasa_vals)),
            "sasa_cv": float(np.std(sasa_vals) / np.mean(sasa_vals))
            if np.mean(sasa_vals) > 0
            else 0,
        }

        line = f"  SASA: {stats['sasa_mean']:.1f} ± {stats['sasa_std']:.1f} A², CV={stats['sasa_cv']:.3f}"
        print(line, file=sys.stderr)
        output_lines.append(line)

        results.append(stats)

    except Exception as e:
        line = f"  ERROR: {e}"
        print(line, file=sys.stderr)
        output_lines.append(line)

# Summary
print("\n" + "=" * 70, file=sys.stderr)
print("SUMMARY", file=sys.stderr)
print("=" * 70, file=sys.stderr)

output_lines.extend(["", "=" * 70, "SUMMARY", "=" * 70])

cham = [r for r in results if r["label"] == "chameleonic"]
non = [r for r in results if r["label"] == "non-chameleonic"]

if cham and non:
    for metric in ["sasa_std", "sasa_range", "sasa_cv"]:
        cham_val = np.mean([r[metric] for r in cham])
        non_val = np.mean([r[metric] for r in non])
        line = f"\n{metric}: Chameleonic={cham_val:.3f}, Non={non_val:.3f}, Gap={cham_val - non_val:+.3f}"
        print(line, file=sys.stderr)
        output_lines.append(line)

# Save
with open("experiments/iter_26_output.txt", "w") as f:
    f.write("\n".join(output_lines))

with open("experiments/iter_26_descriptors.json", "w") as f:
    json.dump(results, f, indent=2)

print("\nResults saved!", file=sys.stderr)
