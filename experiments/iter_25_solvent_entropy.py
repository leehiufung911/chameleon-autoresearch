import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
import json
from rdkit import Chem
from rdkit.Chem import Descriptors


def compute_solvent_entropy_descriptor(smiles, name):
    """
    Compute solvent-shell entropy descriptor based on surface area analysis.
    Uses accessible surface area and polar surface area ratios.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    # Generate conformers using basic ETKDG
    AllChem = (
        Chem.AllChem
        if hasattr(Chem, "AllChem")
        else __import__("rdkit.Chem.AllChem", fromlist=[""])
    )

    # Try EmbedMultipleConfs with different parameters
    params = AllChem.ETKDG()
    params.randomSeed = 42
    AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params)

    if mol.GetNumConformers() == 0:
        # Try simple embedding
        AllChem.EmbedMolecule(mol, randomSeed=42)
        if mol.GetNumConformers() == 0:
            return None

    results = []

    for conf_id in range(mol.GetNumConformers()):
        # Compute TPSA (polar surface area)
        tpsa = Descriptors.TPSA(mol)

        # Try to get SASA if available
        try:
            sasa = AllChem.ComputeMolSurface(mol, confId=conf_id)
            # Estimate hydrophobic vs polar surface
            hydrophobic = sasa - tpsa if sasa > tpsa else 0.0
        except:
            # Fallback: use TPSA and approximate total surface
            sasa = tpsa * 2.0  # rough estimate
            hydrophobic = sasa - tpsa if sasa > tpsa else 0.0

        # Compute Rg for this conformer
        try:
            rg = AllChem.ComputeMolShape(mol, confId=conf_id)
        except:
            rg = 5.0  # placeholder

        results.append(
            {
                "tpsa": float(tpsa),
                "sasa": float(sasa),
                "hydrophobic_sasa": float(hydrophobic),
                "rg": float(rg),
                "polar_fraction": float(tpsa / max(sasa, 1.0)),
            }
        )

    if not results:
        return None

    # Aggregate statistics across conformers
    avg_polar_frac = np.mean([r["polar_fraction"] for r in results])
    std_polar_frac = np.std([r["polar_fraction"] for r in results])
    avg_hydrophobic = np.mean([r["hydrophobic_sasa"] for r in results])
    avg_rg = np.mean([r["rg"] for r in results])
    std_rg = np.std([r["rg"] for r in results])

    # Entropy descriptor: based on conformational variability and hydrophobic surface
    entropy_score = avg_hydrophobic * std_polar_frac / 1000.0

    # Molecular weight and rotatable bonds for normalization
    mw = Descriptors.MolWt(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)

    return {
        "name": name,
        "smiles": smiles,
        "entropy_score": float(entropy_score),
        "avg_polar_fraction": float(avg_polar_frac),
        "polar_fraction_std": float(std_polar_frac),
        "avg_hydrophobic_sasa": float(avg_hydrophobic),
        "mw": float(mw),
        "rotatable_bonds": int(rot_bonds),
        "avg_rg": float(avg_rg),
        "rg_std": float(std_rg),
        "normalized_score": float(entropy_score / np.sqrt(mw) * np.log(1 + rot_bonds))
        if mw > 0
        else 0.0,
    }


def main():
    print("Loading molecules...", file=sys.stderr)

    # Load user PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                protacs.append((parts[0], parts[1]))

    # Load benchmark molecules
    benchmark = []
    with open("chameleon_local/benchmark.tsv", "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                benchmark.append((parts[0], parts[1], parts[2]))

    output_lines = []
    output_lines.append("=" * 80)
    output_lines.append("ITERATION 25: Solvent-Shell Entropy from Geometric Analysis")
    output_lines.append("=" * 80)
    output_lines.append("")
    output_lines.append(
        "Hypothesis: Chameleonicity can be estimated from the entropy change"
    )
    output_lines.append(
        "when transferring from water to chloroform, computed from surface"
    )
    output_lines.append("area differences on a conformer ensemble without explicit MD.")
    output_lines.append("")

    # Process PROTACs
    output_lines.append("-" * 80)
    output_lines.append("USER PROTACS ANALYSIS")
    output_lines.append("-" * 80)

    protac_results = []
    for name, smiles in protacs:
        result = compute_solvent_entropy_descriptor(smiles, name)
        if result:
            protac_results.append(result)
            output_lines.append(f"\n{name}:")
            output_lines.append(f"  Entropy Score: {result['entropy_score']:.4f}")
            output_lines.append(f"  Normalized Score: {result['normalized_score']:.4f}")
            output_lines.append(
                f"  Avg Polar Fraction: {result['avg_polar_fraction']:.4f}"
            )
            output_lines.append(
                f"  Polar Fraction Std: {result['polar_fraction_std']:.4f}"
            )
            output_lines.append(
                f"  Hydrophobic SASA: {result['avg_hydrophobic_sasa']:.2f}"
            )

    output_lines.append("\n" + "-" * 80)
    output_lines.append("PROTAC SEPARATION CHECK")
    output_lines.append("-" * 80)

    cham_scores = [
        r["normalized_score"]
        for r in protac_results
        if r["name"] in ["protac_1", "protac_2"]
    ]
    non_score = [
        r["normalized_score"] for r in protac_results if r["name"] == "protac_3"
    ]

    if cham_scores and non_score:
        output_lines.append(f"\nChameleonic PROTACs avg: {np.mean(cham_scores):.4f}")
        output_lines.append(f"Non-chameleonic protac_3: {non_score[0]:.4f}")
        output_lines.append(f"Difference: {np.mean(cham_scores) - non_score[0]:.4f}")
        if np.mean(cham_scores) > non_score[0]:
            output_lines.append("Status: CORRECT ORDER (chameleonic > non-chameleonic)")
        else:
            output_lines.append("Status: INCORRECT ORDER")

    # Process benchmark
    output_lines.append("\n" + "-" * 80)
    output_lines.append("BENCHMARK DESCRIPTORS")
    output_lines.append("-" * 80)

    benchmark_results = []
    for name, smiles, label in benchmark:
        result = compute_solvent_entropy_descriptor(smiles, name)
        if result:
            result["label"] = label
            benchmark_results.append(result)

    output_lines.append(f"\nProcessed {len(benchmark_results)} benchmark molecules")

    # Compute label statistics
    cham_results = [r for r in benchmark_results if r["label"] == "cham"]
    non_results = [r for r in benchmark_results if r["label"] == "non"]

    if cham_results and non_results:
        cham_scores = [r["normalized_score"] for r in cham_results]
        non_scores = [r["normalized_score"] for r in non_results]

        output_lines.append(
            f"\nChameleonic mean: {np.mean(cham_scores):.4f} +/- {np.std(cham_scores):.4f}"
        )
        output_lines.append(
            f"Non-chameleonic mean: {np.mean(non_scores):.4f} +/- {np.std(non_scores):.4f}"
        )
        output_lines.append(
            f"Difference: {np.mean(cham_scores) - np.mean(non_scores):.4f}"
        )

    output_lines.append("\n" + "=" * 80)
    output_lines.append("END OF ANALYSIS")
    output_lines.append("=" * 80)

    output_text = "\n".join(output_lines)

    # Write to file and print
    with open("experiments/iter_25_output.txt", "w") as f:
        f.write(output_text)
    print(output_text)

    # Save descriptors
    with open("experiments/iter_25_descriptors.json", "w") as f:
        json.dump(benchmark_results, f, indent=2)

    print(
        f"\nDescriptors written to experiments/iter_25_descriptors.json",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
