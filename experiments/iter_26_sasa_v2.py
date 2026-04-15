"""
Iteration 26: SASA Fluctuation Entropy Analysis
"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

print("Starting SASA analysis...", flush=True)

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
import json

print("Imports successful", flush=True)


# Load molecules
def load_molecules(tsv_path):
    """Load molecules from TSV file."""
    print(f"Loading from {tsv_path}...", flush=True)
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"Loaded {len(df)} rows", flush=True)
    molecules = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["SMILES"])
        if mol:
            mol.SetProp("name", str(row["Name"]))
            mol.SetProp("label", str(row.get("Label", "unknown")))
            molecules.append(mol)
    print(f"Parsed {len(molecules)} molecules", flush=True)
    return molecules


def generate_ensembles(mol, n_conformers=20):
    """Generate conformer ensembles."""
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)

    # Minimize
    for conf_idx in range(mol_h.GetNumConformers()):
        AllChem.MMFFOptimizeMolecule(mol_h, confId=conf_idx, mmffVariant="MMFF94s")

    return mol_h


def compute_sasa_ensemble(mol):
    """Compute SASA for all conformers."""
    sasa_values = []
    radii = rdFreeSASA.classifyAtoms(mol)

    for conf_idx in range(mol.GetNumConformers()):
        sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_idx)
        sasa_values.append(sasa)

    return np.array(sasa_values)


def main():
    print("=" * 70, flush=True)
    print("Iteration 26: SASA Fluctuation Entropy", flush=True)
    print("=" * 70, flush=True)

    # Load PROTACs
    print("\nLoading PROTACs...", flush=True)
    protacs = load_molecules("chameleon_local/user_protacs.tsv")

    results = []

    # Process PROTACs
    print("\nProcessing PROTACs:", flush=True)
    print("-" * 70, flush=True)

    for mol in protacs:
        name = mol.GetProp("name")
        label = mol.GetProp("label")
        print(f"\n{name} ({label})", flush=True)

        try:
            # Generate ensemble
            mol_ens = generate_ensembles(mol, n_conformers=20)

            # Compute SASA
            sasa_vals = compute_sasa_ensemble(mol_ens)

            # Stats
            sasa_mean = float(np.mean(sasa_vals))
            sasa_std = float(np.std(sasa_vals))
            sasa_min = float(np.min(sasa_vals))
            sasa_max = float(np.max(sasa_vals))
            sasa_range = sasa_max - sasa_min
            sasa_cv = sasa_std / sasa_mean if sasa_mean > 0 else 0

            print(f"  SASA: {sasa_mean:.1f} ± {sasa_std:.1f} A²", flush=True)
            print(f"  Range: {sasa_range:.1f} A²", flush=True)
            print(f"  CV: {sasa_cv:.3f}", flush=True)

            results.append(
                {
                    "name": name,
                    "label": label,
                    "sasa_mean": sasa_mean,
                    "sasa_std": sasa_std,
                    "sasa_range": sasa_range,
                    "sasa_cv": sasa_cv,
                }
            )
        except Exception as e:
            print(f"  ERROR: {e}", flush=True)

    # Summary
    print("\n" + "=" * 70, flush=True)
    print("Summary", flush=True)
    print("=" * 70, flush=True)

    cham = [r for r in results if r["label"] == "chameleonic"]
    non = [r for r in results if r["label"] == "non-chameleonic"]

    if cham and non:
        for metric in ["sasa_std", "sasa_range", "sasa_cv"]:
            cham_vals = [r[metric] for r in cham]
            non_vals = [r[metric] for r in non]
            print(f"\n{metric}:", flush=True)
            print(f"  Chameleonic: {np.mean(cham_vals):.3f}", flush=True)
            print(f"  Non-chameleonic: {np.mean(non_vals):.3f}", flush=True)
            print(f"  Gap: {np.mean(cham_vals) - np.mean(non_vals):+.3f}", flush=True)

    # Save results
    with open("experiments/iter_26_output.txt", "w") as f:
        for r in results:
            f.write(
                f"{r['name']}: sasa_std={r['sasa_std']:.2f}, sasa_cv={r['sasa_cv']:.3f}\n"
            )

    with open("experiments/iter_26_descriptors.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\nResults saved", flush=True)


if __name__ == "__main__":
    main()
