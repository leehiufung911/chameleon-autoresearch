"""Iteration 26: SASA Fluctuation Entropy Analysis"""

import sys, os
import json
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA

os.environ["PYTHONUNBUFFERED"] = "1"

output = []
output.append("=" * 70)
output.append("ITERATION 26: SASA FLUCTUATION ENTROPY ANALYSIS")
output.append("=" * 70)
output.append("")

# Load PROTACs
df = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t")

results = []

for _, row in df.iterrows():
    name = row["Name"]
    label = row["Label"]
    smiles = row["SMILES"]

    output.append(f"\n{name} ({label})")
    output.append(f"  SMILES: {smiles[:60]}...")

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            output.append("  ERROR: Failed to parse SMILES")
            continue

        mol = Chem.AddHs(mol)

        # Generate small ensemble (10 confs for speed)
        AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=42)

        # Minimize each conformer
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
            "sasa_min": float(np.min(sasa_vals)),
            "sasa_max": float(np.max(sasa_vals)),
            "sasa_range": float(np.max(sasa_vals) - np.min(sasa_vals)),
            "sasa_cv": float(np.std(sasa_vals) / np.mean(sasa_vals))
            if np.mean(sasa_vals) > 0
            else 0,
        }

        output.append(f"  SASA mean: {stats['sasa_mean']:.1f} A²")
        output.append(f"  SASA std:  {stats['sasa_std']:.1f} A²")
        output.append(f"  SASA range: {stats['sasa_range']:.1f} A²")
        output.append(f"  SASA CV: {stats['sasa_cv']:.3f}")

        results.append(stats)

    except Exception as e:
        output.append(f"  ERROR: {e}")

# Summary
output.append("\n" + "=" * 70)
output.append("SUMMARY")
output.append("=" * 70)

cham = [r for r in results if r["label"] == "chameleonic"]
non = [r for r in results if r["label"] == "non-chameleonic"]

if cham and non:
    for metric in ["sasa_std", "sasa_range", "sasa_cv"]:
        cham_val = np.mean([r[metric] for r in cham])
        non_val = np.mean([r[metric] for r in non])
        gap = cham_val - non_val
        output.append(f"\n{metric}:")
        output.append(f"  Chameleonic (protac_1,2): {cham_val:.3f}")
        output.append(f"  Non-chameleonic (protac_3): {non_val:.3f}")
        output.append(f"  Gap: {gap:+.3f}")

# Save to files
output_text = "\n".join(output)

with open("experiments/iter_26_output.txt", "w") as f:
    f.write(output_text)

with open("experiments/iter_26_descriptors.json", "w") as f:
    json.dump(results, f, indent=2)

# Also write to stderr so we can see it
print(output_text, file=sys.stderr)
