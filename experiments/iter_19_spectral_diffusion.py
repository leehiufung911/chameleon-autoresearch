"""
Iteration 19: Graph Spectral Diffusion Analysis

Hypothesis: Chameleonicity can be predicted by analyzing how "efficiently information/energy
flows" through the molecular graph. Molecules with high algebraic connectivity between
distant fragments have folding pathways with low resistance, enabling chameleonic behavior.

Physical Mechanism:
- Model molecule as weighted graph (atoms = nodes, bonds = edges)
- Compute graph Laplacian eigenvalues (spectral decomposition)
- Fiedler value (λ₂) measures "bottleneck-ness" — how well-connected the molecule is
- Spectral entropy captures folding complexity
- Chameleonic molecules have high λ₂ (efficient folding pathways)

This technique is borrowed from network science and spectral graph theory.
"""

import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys

sys.path.insert(0, "chameleon_local")

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import json

print("DEBUG: Script starting...", flush=True)


def build_weighted_graph(mol):
    """
    Build weighted adjacency matrix for the molecular graph.
    Weights account for bond type and atomic properties.
    """
    n_atoms = mol.GetNumAtoms()
    adj = np.zeros((n_atoms, n_atoms))

    # Bond weights: single=1, double=2, triple=3, aromatic=1.5
    bond_weight = {
        Chem.BondType.SINGLE: 1.0,
        Chem.BondType.DOUBLE: 2.0,
        Chem.BondType.TRIPLE: 3.0,
        Chem.BondType.AROMATIC: 1.5,
    }

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        wt = bond_weight.get(bond.GetBondType(), 1.0)

        # Add atomic property modulation
        atom_i = mol.GetAtom(i)
        atom_j = mol.GetAtom(j)

        # Polar atoms (N, O) create stronger "information pathways"
        polar_boost = 1.0
        if atom_i.GetSymbol() in ["N", "O"] or atom_j.GetSymbol() in ["N", "O"]:
            polar_boost = 1.3

        adj[i, j] = wt * polar_boost
        adj[j, i] = wt * polar_boost

    return adj


def compute_graph_laplacian(adj):
    """
    Compute normalized graph Laplacian: L = I - D^(-1/2) @ A @ D^(-1/2)
    """
    # Degree matrix
    degrees = np.sum(adj, axis=1)

    # Avoid division by zero
    degrees = np.where(degrees == 0, 1, degrees)

    # Normalized Laplacian
    d_inv_sqrt = np.diag(1.0 / np.sqrt(degrees))
    laplacian = np.eye(len(adj)) - d_inv_sqrt @ adj @ d_inv_sqrt

    return laplacian


def compute_spectral_descriptors(mol):
    """
    Compute spectral descriptors from graph Laplacian eigenvalues.

    Returns:
        fiedler_value: Second smallest eigenvalue (algebraic connectivity)
        spectral_gap: Difference between largest and smallest eigenvalues
        spectral_entropy: Entropy of eigenvalue distribution
        effective_resistance: Average pairwise resistance (related to commute time)
    """
    adj = build_weighted_graph(mol)
    laplacian = compute_graph_laplacian(adj)

    # Compute eigenvalues
    try:
        eigenvalues = np.linalg.eigvalsh(laplacian)
        eigenvalues = np.sort(eigenvalues)
    except np.linalg.LinAlgError:
        return None

    # Remove zero eigenvalue(s) for connected graphs
    non_zero = eigenvalues[eigenvalues > 1e-10]

    if len(non_zero) < 2:
        return None

    # Fiedler value (algebraic connectivity)
    fiedler = non_zero[0]

    # Spectral gap
    spectral_gap = non_zero[-1] - non_zero[0]

    # Spectral entropy (how "ordered" the spectrum is)
    # Normalize eigenvalues to probabilities
    eig_probs = non_zero / np.sum(non_zero)
    spectral_entropy = -np.sum(eig_probs * np.log(eig_probs + 1e-10))

    # Effective resistance (average commute time between nodes)
    # Related to sum of inverse eigenvalues
    effective_resistance = np.sum(1.0 / non_zero)

    # Size-normalized resistance
    n_nodes = mol.GetNumAtoms()
    normalized_resistance = effective_resistance / (n_nodes**2)

    # Spectral density spread (coefficient of variation)
    spectral_cv = np.std(non_zero) / (np.mean(non_zero) + 1e-10)

    return {
        "fiedler": fiedler,
        "spectral_gap": spectral_gap,
        "spectral_entropy": spectral_entropy,
        "effective_resistance": effective_resistance,
        "normalized_resistance": normalized_resistance,
        "spectral_cv": spectral_cv,
        "n_eigenvalues": len(non_zero),
    }


def compute_diffusion_score(descriptors, mol):
    """
    Combine spectral descriptors into a chameleonicity score.

    Theory: Chameleonic molecules have:
    - Higher Fiedler value (better connected, fewer bottlenecks)
    - Lower spectral entropy (more ordered folding landscape)
    - Lower effective resistance (easier folding pathways)
    """
    if descriptors is None:
        return 0.0, {}

    # Get molecular properties for normalization
    n_atoms = mol.GetNumAtoms()
    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Diffusion efficiency score
    # High Fiedler = good connectivity = easy folding
    # Low entropy = ordered landscape = predictable folding
    # Low resistance = efficient pathways

    fiedler = descriptors["fiedler"]
    entropy = descriptors["spectral_entropy"]
    resistance = descriptors["normalized_resistance"]

    # Combine into diffusion score
    # Normalize by molecule size and flexibility
    if n_rot > 0:
        # Scale by flexibility — flexible molecules need higher scores to be chameleonic
        flexibility_factor = np.log1p(n_rot)

        diffusion_score = (
            fiedler * 2.0  # Connectivity bonus
            + (5.0 - entropy) * 0.5  # Order bonus (entropy typically 2-5)
            + (1.0 / (resistance + 0.01)) * 0.3  # Efficiency bonus
        ) / flexibility_factor
    else:
        diffusion_score = fiedler * 2.0 + (5.0 - entropy) * 0.5

    return diffusion_score, descriptors


def analyze_molecule(smiles, name):
    """Analyze a single molecule"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    desc = compute_spectral_descriptors(mol)
    score, details = compute_diffusion_score(desc, mol)

    return {
        "name": name,
        "smiles": smiles,
        "diffusion_score": score,
        **details,
    }


def main():
    output_lines = []

    output_lines.append("=" * 70)
    output_lines.append("ITERATION 19: Graph Spectral Diffusion Analysis")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append(
        "Hypothesis: Chameleonicity can be predicted from spectral graph"
    )
    output_lines.append("properties that measure 'folding pathway efficiency'.")
    output_lines.append("")
    output_lines.append("Physical Mechanism:")
    output_lines.append("- Fiedler value (λ₂): Measures graph connectivity/bottlenecks")
    output_lines.append("- Spectral entropy: Measures folding landscape order")
    output_lines.append("- Effective resistance: Measures folding pathway efficiency")
    output_lines.append("")
    output_lines.append("This borrows from spectral graph theory in network science.")
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
        result = analyze_molecule(smiles, name)
        if result:
            results.append(result)
            gt = "CHAM" if name in ["protac_1", "protac_2"] else "NON"
            output_lines.append(f"{name} ({gt}):")
            output_lines.append(f"  Fiedler (λ₂): {result.get('fiedler', 0):.4f}")
            output_lines.append(
                f"  Spectral entropy: {result.get('spectral_entropy', 0):.3f}"
            )
            output_lines.append(
                f"  Norm. resistance: {result.get('normalized_resistance', 0):.4f}"
            )
            output_lines.append(f"  Spectral CV: {result.get('spectral_cv', 0):.3f}")
            output_lines.append(f"  Diffusion Score: {result['diffusion_score']:.4f}")
            output_lines.append("")

    # Check separation
    cham_scores = [
        r["diffusion_score"] for r in results if r["name"] in ["protac_1", "protac_2"]
    ]
    non_scores = [r["diffusion_score"] for r in results if r["name"] == "protac_3"]

    if cham_scores and non_scores:
        output_lines.append("-" * 70)
        output_lines.append("SEPARATION ANALYSIS")
        output_lines.append("-" * 70)
        output_lines.append("")
        output_lines.append(
            f"Chameleonic PROTACs avg score: {np.mean(cham_scores):.4f}"
        )
        output_lines.append(f"Non-chameleonic avg score: {np.mean(non_scores):.4f}")
        output_lines.append(
            f"Difference: {np.mean(cham_scores) - np.mean(non_scores):+.4f}"
        )
        output_lines.append("")

        if np.mean(cham_scores) > np.mean(non_scores):
            output_lines.append("Result: CORRECT ORDERING")
            output_lines.append(
                "Diffusion score successfully separates chameleonic from non-chameleonic"
            )
        else:
            output_lines.append("Result: FAILED")
            output_lines.append("Diffusion score does not separate correctly")
        output_lines.append("")

    # Analyze benchmark
    output_lines.append("-" * 70)
    output_lines.append("BENCHMARK ANALYSIS (49 molecules)")
    output_lines.append("-" * 70)
    output_lines.append("")

    benchmark = []
    with open("chameleon_local/benchmark.tsv") as f:
        header = next(f).strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                name = parts[0]
                smiles = parts[1]
                label = parts[2]
                benchmark.append((name, smiles, label))

    benchmark_results = []
    for name, smiles, label in benchmark:
        result = analyze_molecule(smiles, name)
        if result:
            result["label"] = label
            benchmark_results.append(result)

    # Compute statistics
    cham_values = [
        r["diffusion_score"] for r in benchmark_results if r["label"] == "chameleonic"
    ]
    noncham_values = [
        r["diffusion_score"]
        for r in benchmark_results
        if r["label"] == "non_chameleonic"
    ]

    if cham_values and noncham_values:
        output_lines.append(
            f"Chameleonic (n={len(cham_values)}): mean={np.mean(cham_values):.4f}, std={np.std(cham_values):.4f}"
        )
        output_lines.append(
            f"Non-chameleonic (n={len(noncham_values)}): mean={np.mean(noncham_values):.4f}, std={np.std(noncham_values):.4f}"
        )
        output_lines.append(
            f"Separation: {np.mean(cham_values) - np.mean(noncham_values):+.4f}"
        )
        output_lines.append("")

    # Show sample molecules
    output_lines.append("Sample chameleonic molecules:")
    for r in sorted(
        [br for br in benchmark_results if br["label"] == "chameleonic"],
        key=lambda x: x["diffusion_score"],
        reverse=True,
    )[:5]:
        output_lines.append(f"  {r['name']}: {r['diffusion_score']:.4f}")
    output_lines.append("")

    output_lines.append("Sample non-chameleonic molecules:")
    for r in sorted(
        [br for br in benchmark_results if br["label"] == "non_chameleonic"],
        key=lambda x: x["diffusion_score"],
        reverse=True,
    )[:5]:
        output_lines.append(f"  {r['name']}: {r['diffusion_score']:.4f}")
    output_lines.append("")

    output_lines.append("=" * 70)
    output_lines.append("PHYSICAL INTERPRETATION")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append(
        "Graph Spectral Diffusion Analysis treats the molecule as a network"
    )
    output_lines.append("where information/energy flows along bonds.")
    output_lines.append("")
    output_lines.append("Key insights:")
    output_lines.append(
        "1. Fiedler value (λ₂): Higher = better connectivity = easier folding"
    )
    output_lines.append(
        "2. Spectral entropy: Lower = more ordered landscape = predictable folding"
    )
    output_lines.append(
        "3. Effective resistance: Lower = efficient pathways = fast folding"
    )
    output_lines.append("")
    output_lines.append(
        "PEG linkers create multiple parallel pathways (high connectivity),"
    )
    output_lines.append("while alkyl linkers create bottlenecks (low connectivity).")
    output_lines.append("")
    output_lines.append(
        "This is a purely 2D topological approach from network science,"
    )
    output_lines.append("avoiding all 3D conformer bias issues.")

    output = "\n".join(output_lines)

    # Write output
    with open("experiments/iter_19_output.txt", "w") as f:
        f.write(output)

    print(output, flush=True)

    # Write descriptors
    with open("experiments/iter_19_descriptors.json", "w") as f:
        json.dump(benchmark_results, f, indent=2)

    print(
        f"\nWrote {len(benchmark_results)} molecule descriptors to iter_19_descriptors.json",
        flush=True,
    )


if __name__ == "__main__":
    main()
