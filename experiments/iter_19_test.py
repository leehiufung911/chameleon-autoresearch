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


def build_weighted_graph(mol):
    """Build weighted adjacency matrix for the molecular graph."""
    n_atoms = mol.GetNumAtoms()
    adj = np.zeros((n_atoms, n_atoms))

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

        atom_i = mol.GetAtom(i)
        atom_j = mol.GetAtom(j)

        polar_boost = 1.0
        if atom_i.GetSymbol() in ["N", "O"] or atom_j.GetSymbol() in ["N", "O"]:
            polar_boost = 1.3

        adj[i, j] = wt * polar_boost
        adj[j, i] = wt * polar_boost

    return adj


def compute_graph_laplacian(adj):
    """Compute normalized graph Laplacian."""
    degrees = np.sum(adj, axis=1)
    degrees = np.where(degrees == 0, 1, degrees)

    d_inv_sqrt = np.diag(1.0 / np.sqrt(degrees))
    laplacian = np.eye(len(adj)) - d_inv_sqrt @ adj @ d_inv_sqrt

    return laplacian


def compute_spectral_descriptors(mol):
    """Compute spectral descriptors from graph Laplacian eigenvalues."""
    adj = build_weighted_graph(mol)
    laplacian = compute_graph_laplacian(adj)

    try:
        eigenvalues = np.linalg.eigvalsh(laplacian)
        eigenvalues = np.sort(eigenvalues)
    except np.linalg.LinAlgError:
        return None

    non_zero = eigenvalues[eigenvalues > 1e-10]

    if len(non_zero) < 2:
        return None

    fiedler = non_zero[0]
    spectral_gap = non_zero[-1] - non_zero[0]

    eig_probs = non_zero / np.sum(non_zero)
    spectral_entropy = -np.sum(eig_probs * np.log(eig_probs + 1e-10))

    effective_resistance = np.sum(1.0 / non_zero)
    n_nodes = mol.GetNumAtoms()
    normalized_resistance = effective_resistance / (n_nodes**2)

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
    """Combine spectral descriptors into a chameleonicity score."""
    if descriptors is None:
        return 0.0, {}

    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)

    fiedler = descriptors["fiedler"]
    entropy = descriptors["spectral_entropy"]
    resistance = descriptors["normalized_resistance"]

    if n_rot > 0:
        flexibility_factor = np.log1p(n_rot)
        diffusion_score = (
            fiedler * 2.0 + (5.0 - entropy) * 0.5 + (1.0 / (resistance + 0.01)) * 0.3
        ) / flexibility_factor
    else:
        diffusion_score = fiedler * 2.0 + (5.0 - entropy) * 0.5

    return diffusion_score, descriptors


def analyze_molecule(smiles, name):
    """Analyze a single molecule."""
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
    output_lines.append("Hypothesis: Chameleonicity from spectral graph properties")
    output_lines.append("that measure 'folding pathway efficiency'.")
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
            output_lines.append(f"  Fiedler: {result.get('fiedler', 0):.4f}")
            output_lines.append(f"  Entropy: {result.get('spectral_entropy', 0):.3f}")
            output_lines.append(
                f"  Resistance: {result.get('normalized_resistance', 0):.4f}"
            )
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
        output_lines.append(f"Chameleonic PROTACs avg: {np.mean(cham_scores):.4f}")
        output_lines.append(f"Non-chameleonic avg: {np.mean(non_scores):.4f}")
        output_lines.append(
            f"Difference: {np.mean(cham_scores) - np.mean(non_scores):+.4f}"
        )
        output_lines.append("")

        if np.mean(cham_scores) > np.mean(non_scores):
            output_lines.append("Result: CORRECT ORDERING - SUCCESS")
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

    output_lines.append("=" * 70)
    output_lines.append("PHYSICAL INTERPRETATION")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append("Graph Spectral Analysis from network science:")
    output_lines.append(
        "1. Fiedler value: Higher = better connectivity = easier folding"
    )
    output_lines.append(
        "2. Spectral entropy: Lower = more ordered = predictable folding"
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
        "Verdict: This is a purely 2D topological approach from network"
    )
    output_lines.append("science, avoiding all 3D conformer bias issues.")

    output = "\n".join(output_lines)

    # Write output
    with open("experiments/iter_19_output.txt", "w") as f:
        f.write(output)

    print(output, flush=True)

    # Write descriptors
    with open("experiments/iter_19_descriptors.json", "w") as f:
        json.dump(benchmark_results, f, indent=2)

    print(
        f"\nWrote {len(benchmark_results)} descriptors to iter_19_descriptors.json",
        flush=True,
    )


if __name__ == "__main__":
    main()
