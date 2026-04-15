"""
Iteration 26: SASA Fluctuation Entropy Analysis - NEW VERSION
"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
import json

print("SASA ENTROPY ANALYSIS - Starting", flush=True)


# Load molecules
def load_molecules(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t")
    molecules = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["SMILES"])
        if mol:
            mol.SetProp("name", str(row["Name"]))
            mol.SetProp("label", str(row.get("Label", "unknown")))
            molecules.append(mol)
    return molecules


def generate_ensemble(mol, n_conformers=20):
    """Generate minimized conformer ensemble."""
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(
        mol_h, numConfs=n_conformers, randomSeed=42, params=AllChem.ETKDGv3()
    )
    for i in range(mol_h.GetNumConformers()):
        AllChem.MMFFOptimizeMolecule(mol_h, confId=i, mmffVariant="MMFF94s")
    return mol_h


def compute_sasa_stats(mol):
    """Compute SASA statistics across ensemble."""
    radii = rdFreeSASA.classifyAtoms(mol)
    sasa_vals = []
    for i in range(mol.GetNumConformers()):
        sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=i)
        sasa_vals.append(sasa)
    sasa_vals = np.array(sasa_vals)

    return {
        "sasa_mean": float(np.mean(sasa_vals)),
        "sasa_std": float(np.std(sasa_vals)),
        "sasa_min": float(np.min(sasa_vals)),
        "sasa_max": float(np.max(sasa_vals)),
        "sasa_range": float(np.max(sasa_vals) - np.min(sasa_vals)),
        "sasa_cv": float(np.std(sasa_vals) / np.mean(sasa_vals))
        if np.mean(sasa_vals) > 0
        else 0,
    }


def main():
    # Load PROTACs
    protacs = load_molecules("chameleon_local/user_protacs.tsv")

    print("\n" + "=" * 70)
    print("ITERATION 26: SASA FLUCTUATION ENTROPY")
    print("=" * 70)

    results = []

    for mol in protacs:
        name = mol.GetProp("name")
        label = mol.GetProp("label")

        print(f"\n{name} ({label})", flush=True)

        try:
            mol_ens = generate_ensemble(mol, n_conformers=25)
            stats = compute_sasa_stats(mol_ens)

            print(f"  SASA mean: {stats['sasa_mean']:.1f} A^2", flush=True)
            print(f"  SASA std:  {stats['sasa_std']:.1f} A^2", flush=True)
            print(f"  SASA range: {stats['sasa_range']:.1f} A^2", flush=True)
            print(f"  SASA CV: {stats['sasa_cv']:.3f}", flush=True)

            results.append({"name": name, "label": label, **stats})
        except Exception as e:
            print(f"  ERROR: {e}", flush=True)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    cham = [r for r in results if r["label"] == "chameleonic"]
    non = [r for r in results if r["label"] == "non-chameleonic"]

    for metric in ["sasa_std", "sasa_range", "sasa_cv"]:
        if cham and non:
            cham_val = np.mean([r[metric] for r in cham])
            non_val = np.mean([r[metric] for r in non])
            print(f"\n{metric}:")
            print(f"  Chameleonic:    {cham_val:.3f}")
            print(f"  Non-chameleonic: {non_val:.3f}")
            print(f"  Gap: {cham_val - non_val:+.3f}")

    # Save
    with open("experiments/iter_26_output.txt", "w") as f:
        for r in results:
            f.write(
                f"{r['name']}: sasa_std={r['sasa_std']:.2f}, sasa_cv={r['sasa_cv']:.3f}\n"
            )

    with open("experiments/iter_26_descriptors.json", "w") as f:
        json.dump(results, f, indent=2)

    print(
        "\nResults saved to experiments/iter_26_output.txt and iter_26_descriptors.json"
    )


if __name__ == "__main__":
    main()
