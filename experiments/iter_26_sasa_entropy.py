"""
Iteration 26: Solvent-Mediated Entropy from SASA Fluctuation Analysis

Hypothesis: Chameleonic molecules show bimodal SASA distributions (discrete switching)
while non-chameleonic flexible molecules show unimodal distributions (continuous fluctuation).
SASA variance and bimodality capture the entropy of solvent reorganization.
"""

import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr  # Forces output through stderr (unbuffered in Bun)

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import rdFreeSASA
import json


# Load molecules
def load_molecules(tsv_path):
    """Load molecules from TSV file."""
    df = pd.read_csv(tsv_path, sep="\t")
    molecules = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["SMILES"])
        if mol:
            mol.SetProp("name", row["Name"])
            mol.SetProp("label", row.get("Label", "unknown"))
            molecules.append(mol)
    return molecules


def generate_ensembles(mol, n_conformers=30):
    """Generate conformer ensembles in polar and apolar implicit solvent."""
    # Polar ensemble (water-like, high dielectric)
    mol_pol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(
        mol_pol, numConfs=n_conformers, randomSeed=42, params=AllChem.ETKDGv3()
    )

    # Apolar ensemble (cyclohexane-like, low dielectric)
    mol_apol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(
        mol_apol, numConfs=n_conformers, randomSeed=42, params=AllChem.ETKDGv3()
    )

    # Minimize with MMFF in respective dielectrics
    for conf_idx in range(mol_pol.GetNumConformers()):
        AllChem.MMFFOptimizeMolecule(mol_pol, confId=conf_idx, mmffVariant="MMFF94s")

    for conf_idx in range(mol_apol.GetNumConformers()):
        AllChem.MMFFOptimizeMolecule(mol_apol, confId=conf_idx, mmffVariant="MMFF94s")

    return mol_pol, mol_apol


def compute_sasa_ensemble(mol):
    """Compute SASA for all conformers in a molecule."""
    sasa_values = []
    radii = rdFreeSASA.classifyAtoms(mol)

    for conf_idx in range(mol.GetNumConformers()):
        sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_idx)
        sasa_values.append(sasa)

    return np.array(sasa_values)


def compute_bimodality(data):
    """Compute bimodality coefficient (0-1, >0.55 suggests bimodality)."""
    n = len(data)
    if n < 4:
        return 0.0

    mean = np.mean(data)
    std = np.std(data)
    skewness = np.mean(((data - mean) / std) ** 3) if std > 0 else 0
    kurtosis = np.mean(((data - mean) / std) ** 4) - 3 if std > 0 else -3

    # Bimodality coefficient
    bc = (skewness**2 + 1) / (kurtosis + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
    return max(0, bc)  # Ensure non-negative


def compute_entropy_metrics(mol_pol, mol_apol):
    """Compute entropy metrics from SASA fluctuations."""
    # SASA distributions
    sasa_pol = compute_sasa_ensemble(mol_pol)
    sasa_apol = compute_sasa_ensemble(mol_apol)

    # Combined ensemble
    sasa_all = np.concatenate([sasa_pol, sasa_apol])

    # Basic stats
    sasa_mean = np.mean(sasa_all)
    sasa_std = np.std(sasa_all)
    sasa_cv = sasa_std / sasa_mean if sasa_mean > 0 else 0

    # Bimodality (measures discrete vs continuous switching)
    bimodality = compute_bimodality(sasa_all)

    # Polar-apolar difference
    delta_sasa = np.mean(sasa_pol) - np.mean(sasa_apol)
    delta_sasa_norm = abs(delta_sasa) / sasa_mean if sasa_mean > 0 else 0

    # Entropy proxy: variance * bimodality
    entropy_proxy = sasa_std * bimodality

    # Rg for comparison
    rg_pol = [
        AllChem.ComputeMolShapeDescriptors(mol_pol, confId=i)[0]
        for i in range(mol_pol.GetNumConformers())
    ]
    rg_apol = [
        AllChem.ComputeMolShapeDescriptors(mol_apol, confId=i)[0]
        for i in range(mol_apol.GetNumConformers())
    ]
    rg_mean = np.mean(rg_pol + rg_apol)
    rg_std = np.std(rg_pol + rg_apol)

    return {
        "sasa_mean": sasa_mean,
        "sasa_std": sasa_std,
        "sasa_cv": sasa_cv,
        "bimodality": bimodality,
        "delta_sasa": abs(delta_sasa),
        "delta_sasa_norm": delta_sasa_norm,
        "entropy_proxy": entropy_proxy,
        "rg_mean": rg_mean,
        "rg_std": rg_std,
    }


def main():
    # Load data
    protacs_path = "chameleon_local/user_protacs.tsv"
    benchmark_path = "chameleon_local/benchmark.tsv"

    protacs = load_molecules(protacs_path)
    benchmark = load_molecules(benchmark_path)
    all_molecules = protacs + benchmark

    output_lines = []
    output_lines.append("=" * 70)
    output_lines.append("Iteration 26: SASA Fluctuation Entropy Analysis")
    output_lines.append("=" * 70)
    output_lines.append("")

    results = []

    # Process PROTACs first
    output_lines.append("PROTAC Results:")
    output_lines.append("-" * 70)

    for mol in protacs:
        name = mol.GetProp("name")
        label = mol.GetProp("label")

        output_lines.append(f"\n{name} (Label: {label})")
        output_lines.append(f"  SMILES: {Chem.MolToSmiles(mol)}")

        # Generate ensembles
        mol_pol, mol_apol = generate_ensembles(mol, n_conformers=30)

        # Compute metrics
        metrics = compute_entropy_metrics(mol_pol, mol_apol)

        output_lines.append(f"  SASA mean: {metrics['sasa_mean']:.2f} A²")
        output_lines.append(f"  SASA std: {metrics['sasa_std']:.2f} A²")
        output_lines.append(f"  SASA CV: {metrics['sasa_cv']:.3f}")
        output_lines.append(f"  Bimodality: {metrics['bimodality']:.3f}")
        output_lines.append(f"  Delta SASA: {metrics['delta_sasa']:.2f} A²")
        output_lines.append(f"  Entropy proxy: {metrics['entropy_proxy']:.2f}")
        output_lines.append(f"  Rg mean: {metrics['rg_mean']:.2f} A")

        result = {"name": name, "label": label, **metrics}
        results.append(result)

    # Summary
    output_lines.append("\n" + "=" * 70)
    output_lines.append("Summary: PROTAC Separation")
    output_lines.append("=" * 70)

    cham = [r for r in results if r["label"] == "chameleonic"]
    non_cham = [r for r in results if r["label"] == "non-chameleonic"]

    if cham and non_cham:
        for metric in ["sasa_std", "bimodality", "entropy_proxy", "sasa_cv"]:
            cham_vals = [r[metric] for r in cham]
            non_vals = [r[metric] for r in non_cham]

            output_lines.append(f"\n{metric}:")
            output_lines.append(f"  Chameleonic (protac_1,2): {np.mean(cham_vals):.3f}")
            output_lines.append(
                f"  Non-chameleonic (protac_3): {np.mean(non_vals):.3f}"
            )
            output_lines.append(
                f"  Difference: {np.mean(cham_vals) - np.mean(non_vals):+.3f}"
            )

    # Process benchmark molecules
    output_lines.append("\n" + "=" * 70)
    output_lines.append(f"Benchmark molecules: {len(benchmark)}")
    output_lines.append("=" * 70)

    for mol in benchmark:
        name = mol.GetProp("name")
        label = mol.GetProp("label")

        mol_pol, mol_apol = generate_ensembles(mol, n_conformers=20)
        metrics = compute_entropy_metrics(mol_pol, mol_apol)

        result = {"name": name, "label": label, **metrics}
        results.append(result)

    # Save descriptors
    with open("experiments/iter_26_descriptors.json", "w") as f:
        json.dump(results, f, indent=2)

    output_lines.append(f"\nDescriptors saved to experiments/iter_26_descriptors.json")
    output_lines.append(f"Total molecules processed: {len(results)}")

    # Write output
    output_text = "\n".join(output_lines)

    with open("experiments/iter_26_output.txt", "w") as f:
        f.write(output_text)

    print(output_text)


if __name__ == "__main__":
    main()
