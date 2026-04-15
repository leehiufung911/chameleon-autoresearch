"""Iteration 26: SASA Fluctuation Entropy Analysis - FINAL VERSION"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

print("=" * 70)
print("ITERATION 26: SASA FLUCTUATION ENTROPY ANALYSIS")
print("=" * 70)

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
import json

print("\nLoading PROTACs...")

# Load PROTACs
df = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t")

results = []

for idx, row in df.iterrows():
    name = row["Name"]
    label = row["Label"]
    smiles = row["SMILES"]

    print(f"\n{name} ({label})")
    print(f"  SMILES: {smiles[:60]}...")

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print("  ERROR: Failed to parse")
            continue

        mol = Chem.AddHs(mol)

        # Generate ensemble
        print("  Generating conformers...")
        AllChem.EmbedMultipleConfs(mol, numConfs=12, randomSeed=42)

        # Minimize
        print("  Minimizing...")
        for i in range(mol.GetNumConformers()):
            AllChem.MMFFOptimizeMolecule(mol, confId=i)

        # Compute SASA
        print("  Computing SASA...")
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

        print(f"  SASA: {stats['sasa_mean']:.1f} ± {stats['sasa_std']:.1f} A²")
        print(f"  Range: {stats['sasa_range']:.1f} A²")
        print(f"  CV: {stats['sasa_cv']:.3f}")

        results.append(stats)

    except Exception as e:
        print(f"  ERROR: {e}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

cham = [r for r in results if r["label"] == "chameleonic"]
non = [r for r in results if r["label"] == "non-chameleonic"]

if cham and non:
    for metric in ["sasa_std", "sasa_range", "sasa_cv"]:
        cham_val = np.mean([r[metric] for r in cham])
        non_val = np.mean([r[metric] for r in non])
        gap = cham_val - non_val
        print(f"\n{metric}:")
        print(f"  Chameleonic: {cham_val:.3f}")
        print(f"  Non-chameleonic: {non_val:.3f}")
        print(f"  Gap: {gap:+.3f}")

# Save results
print("\n" + "=" * 70)
print("Saving results...")

output_lines = []
for r in results:
    output_lines.append(
        f"{r['name']}: sasa_std={r['sasa_std']:.2f}, sasa_cv={r['sasa_cv']:.3f}"
    )

with open("experiments/iter_26_output.txt", "w") as f:
    f.write("\n".join(output_lines))

with open("experiments/iter_26_descriptors.json", "w") as f:
    json.dump(results, f, indent=2)

print("Results saved to experiments/iter_26_output.txt")
print("Descriptors saved to experiments/iter_26_descriptors.json")
