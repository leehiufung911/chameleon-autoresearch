"""
Iteration 20: Contact Number Entropy Analysis

Hypothesis: Chameleonic molecules have bimodal contact number distributions
due to two-state behavior (collapsed vs extended), while non-chameleonic
molecules have unimodal distributions.

Physical Mechanism (from polymer physics):
- Contact number (Nc): count of atom pairs within topological distance cutoff
- Chameleonic: Two populations -> high entropy of Nc distribution
- Non-chameleonic: Single population -> low entropy of Nc distribution

This is a 2D topological approach measuring "topological folding diversity"
rather than 3D conformer sampling.
"""

import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys
import json
import math
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def compute_topological_contact_entropy(mol, cutoff=8):
    """
    Compute topological contact number entropy.

    Instead of 3D distances, use graph topological distances (shortest path length).
    This measures "folding pathway diversity" in the molecular graph topology.

    Args:
        mol: RDKit molecule
        cutoff: maximum topological distance to count as "contact"

    Returns:
        dict with contact entropy metrics
    """
    if mol is None:
        return None

    n_atoms = mol.GetNumAtoms()
    if n_atoms < 2:
        return None

    # Compute shortest path distances (topological distances)
    # Use Floyd-Warshall or repeated BFS
    distances = {}

    for start_atom in range(n_atoms):
        # BFS from start_atom
        visited = {start_atom: 0}
        queue = [start_atom]

        while queue:
            current = queue.pop(0)
            current_dist = visited[current]

            atom = mol.GetAtomWithIdx(current)
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in visited:
                    visited[n_idx] = current_dist + 1
                    queue.append(n_idx)

        distances[start_atom] = visited

    # Build "contact" distribution by sampling random pairs
    # In polymer physics, contacts are pairs within certain distance
    contact_numbers = []

    # Sample multiple "folding configurations" by varying focus
    # Each configuration counts contacts from a different starting point
    n_samples = min(20, n_atoms)

    for start_idx in range(0, n_atoms, max(1, n_atoms // n_samples)):
        n_contacts = 0

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                # Skip bonded pairs (distance = 1)
                if distances[i][j] > 1 and distances[i][j] <= cutoff:
                    n_contacts += 1

        contact_numbers.append(n_contacts)

    if len(contact_numbers) < 2:
        return None

    # Compute statistics
    nc_mean = sum(contact_numbers) / len(contact_numbers)
    nc_min = min(contact_numbers)
    nc_max = max(contact_numbers)
    nc_range = nc_max - nc_min

    # Standard deviation
    variance = sum((x - nc_mean) ** 2 for x in contact_numbers) / len(contact_numbers)
    nc_std = math.sqrt(variance)

    # Shannon entropy (binning approach)
    n_bins = 5
    if nc_range > 0:
        bin_width = nc_range / n_bins
        bins = [0] * n_bins
        for x in contact_numbers:
            bin_idx = min(int((x - nc_min) / bin_width), n_bins - 1)
            bins[bin_idx] += 1
    else:
        bins = [len(contact_numbers)] + [0] * (n_bins - 1)

    total = sum(bins)
    if total > 0:
        entropy = 0.0
        for count in bins:
            if count > 0:
                p = count / total
                entropy -= p * math.log(p)
    else:
        entropy = 0.0

    # Coefficient of variation
    cv = nc_std / nc_mean if nc_mean > 0 else 0.0

    # Bimodality coefficient
    if nc_std > 0 and len(contact_numbers) > 3:
        skewness = (
            sum((x - nc_mean) ** 3 for x in contact_numbers)
            / len(contact_numbers)
            / (nc_std**3)
        )
        kurtosis = (
            sum((x - nc_mean) ** 4 for x in contact_numbers)
            / len(contact_numbers)
            / (nc_std**4)
            - 3
        )
        bimodality = (skewness**2 + 1) / (kurtosis + 3) if (kurtosis + 3) != 0 else 0
    else:
        skewness = kurtosis = bimodality = 0.0

    return {
        "nc_mean": nc_mean,
        "nc_range": nc_range,
        "nc_std": nc_std,
        "entropy": entropy,
        "cv": cv,
        "bimodality": bimodality,
        "skewness": skewness,
    }


def main():
    output_lines = []

    output_lines.append("=" * 70)
    output_lines.append("ITERATION 20: Contact Number Entropy Analysis")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append("Hypothesis: Chameleonic molecules have high topological")
    output_lines.append("contact number entropy (bimodal distributions) while")
    output_lines.append("non-chameleonic molecules have low entropy (unimodal).")
    output_lines.append("")

    # Read user PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                protacs.append((parts[0], parts[1]))

    results = []

    output_lines.append("-" * 70)
    output_lines.append("USER PROTACS ANALYSIS")
    output_lines.append("-" * 70)
    output_lines.append("")

    for name, smiles in protacs:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol = Chem.AddHs(mol)

        result = compute_topological_contact_entropy(mol, cutoff=8)
        if result:
            result["name"] = name
            result["smiles"] = smiles
            results.append(result)

            gt = "CHAM" if name in ["protac_1", "protac_2"] else "NON"
            output_lines.append(f"{name} ({gt}):")
            output_lines.append(f"  Contact mean: {result['nc_mean']:.1f}")
            output_lines.append(f"  Contact range: {result['nc_range']:.1f}")
            output_lines.append(f"  Entropy: {result['entropy']:.4f}")
            output_lines.append(f"  CV: {result['cv']:.4f}")
            output_lines.append(f"  Bimodality: {result['bimodality']:.4f}")
            output_lines.append("")

    # Check separation
    cham_entropy = [
        r["entropy"] for r in results if r["name"] in ["protac_1", "protac_2"]
    ]
    non_entropy = [r["entropy"] for r in results if r["name"] == "protac_3"]

    if cham_entropy and non_entropy:
        output_lines.append("-" * 70)
        output_lines.append("SEPARATION ANALYSIS")
        output_lines.append("-" * 70)
        output_lines.append("")

        avg_cham = sum(cham_entropy) / len(cham_entropy)
        avg_non = sum(non_entropy) / len(non_entropy)

        output_lines.append(f"Chameleonic avg entropy: {avg_cham:.4f}")
        output_lines.append(f"Non-chameleonic entropy: {avg_non:.4f}")
        output_lines.append(f"Difference: {avg_cham - avg_non:+.4f}")
        output_lines.append("")

        if avg_cham > avg_non:
            output_lines.append(
                "Result: CORRECT ORDERING - Chameleonic > Non-chameleonic"
            )
        else:
            output_lines.append("Result: FAILED - Wrong ordering")
        output_lines.append("")

    # Analyze benchmark
    output_lines.append("-" * 70)
    output_lines.append("BENCHMARK ANALYSIS (49 molecules)")
    output_lines.append("-" * 70)
    output_lines.append("")

    benchmark = []
    with open("chameleon_local/benchmark.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                benchmark.append((parts[0], parts[1]))

    output_lines.append(f"Loaded {len(benchmark)} benchmark molecules")
    output_lines.append("")

    benchmark_results = []
    for name, smiles in benchmark:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol = Chem.AddHs(mol)

        result = compute_topological_contact_entropy(mol, cutoff=8)
        if result:
            result["name"] = name
            benchmark_results.append(result)
            output_lines.append(
                f"{name:25s} Entropy={result['entropy']:.4f}  CV={result['cv']:.4f}"
            )

    output_lines.append("")
    output_lines.append(f"Successfully analyzed {len(benchmark_results)} molecules")
    output_lines.append("")

    # Write JSON descriptors
    with open("experiments/iter_20_descriptors.json", "w") as f:
        json.dump(benchmark_results, f, indent=2)

    output_lines.append("=" * 70)
    output_lines.append("PHYSICAL INTERPRETATION")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append("Contact Number Entropy from Polymer Physics:")
    output_lines.append("- Nc: Count of atom pairs within topological distance cutoff")
    output_lines.append("- Entropy: Shannon entropy of Nc distribution")
    output_lines.append("- High entropy = diverse folding configurations (chameleonic)")
    output_lines.append(
        "- Low entropy = single dominant configuration (non-chameleonic)"
    )
    output_lines.append("")
    output_lines.append("This is a 2D topological approach from polymer physics.")
    output_lines.append("Measures 'topological folding diversity' without 3D confs.")

    output = "\n".join(output_lines)

    # Write output
    with open("experiments/iter_20_output.txt", "w") as f:
        f.write(output)

    print(
        f"Analyzed {len(results)} PROTACs and {len(benchmark_results)} benchmark molecules"
    )
    print(f"Output written to experiments/iter_20_output.txt")
    print(f"Descriptors written to experiments/iter_20_descriptors.json")


if __name__ == "__main__":
    main()
