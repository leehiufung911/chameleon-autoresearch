#!/usr/bin/env python3
"""
Iteration 24: Molecular Spectral Flexibility Index
A novel approach from solid-state physics and protein dynamics
"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import json


def compute_graph_spectral_properties(smiles):
    """Compute spectral properties from molecular graph."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    n_atoms = mol.GetNumAtoms()

    # Build adjacency matrix
    adj = np.zeros((n_atoms, n_atoms))
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if bond.GetIsAromatic():
            w = 1.5
        else:
            bt = bond.GetBondType()
            w = (
                1.0
                if bt == Chem.BondType.SINGLE
                else 2.0
                if bt == Chem.BondType.DOUBLE
                else 3.0
            )
        adj[i, j] = w
        adj[j, i] = w

    # Laplacian: L = D - A
    degrees = np.sum(adj, axis=1)
    laplacian = np.diag(degrees) - adj

    # Compute eigenvalues
    try:
        eigvals = np.linalg.eigvalsh(laplacian)
        nonzero = eigvals[eigvals > 1e-10]
        if len(nonzero) == 0:
            return None

        spectral_width = np.max(nonzero) / np.min(nonzero)
        spectral_gap = nonzero[1] - nonzero[0] if len(nonzero) > 1 else 0
        algebraic_connectivity = nonzero[0]

    except:
        return None

    # Additional descriptors
    n_rot = Descriptors.NumRotatableBonds(mol)
    mw = Descriptors.MolWt(mol)

    # Normalized metrics
    flex_index = spectral_width / (n_rot + 1)  # Flexibility per rotatable bond

    return {
        "spectral_width": float(spectral_width),
        "spectral_gap": float(spectral_gap),
        "algebraic_connectivity": float(algebraic_connectivity),
        "flexibility_index": float(flex_index),
        "n_atoms": n_atoms,
        "n_rot_bonds": n_rot,
        "mw": float(mw),
    }


def main():
    # Initialize output
    lines = []

    lines.append("=" * 70)
    lines.append("ITERATION 24: Molecular Spectral Flexibility Index")
    lines.append("Graph spectral analysis from solid-state physics")
    lines.append("=" * 70)
    lines.append("")

    # User PROTACs
    user_protacs = [
        (
            "protac_1",
            "O=C(C(N1C(C2=CC=CC(NCCOCCOCCOCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=C2C1=O)=O)CC7)NC7=O",
        ),
        (
            "protac_2",
            "O=C(C(N1C(C2=CC=CC(NC(COCCOCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
        ),
        (
            "protac_3",
            "O=C(C(N1C(C2=CC=CC(OCC(NCCCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
        ),
    ]

    # Analyze user PROTACs
    lines.append("USER PROTAC ANALYSIS")
    lines.append("-" * 50)

    user_results = []
    for name, smiles in user_protacs:
        props = compute_graph_spectral_properties(smiles)
        if props:
            user_results.append({"name": name, **props})
            lines.append(f"\n{name}:")
            lines.append(f"  Spectral width: {props['spectral_width']:.4f}")
            lines.append(f"  Spectral gap: {props['spectral_gap']:.4f}")
            lines.append(f"  Flexibility index: {props['flexibility_index']:.6f}")
            lines.append(f"  Rotatable bonds: {props['n_rot_bonds']}")

    # Compare
    if len(user_results) == 3:
        lines.append("")
        lines.append("=" * 50)
        lines.append("PROTAC COMPARISON")
        lines.append("=" * 50)

        cham = [r for r in user_results if r["name"] in ["protac_1", "protac_2"]]
        non_cham = [r for r in user_results if r["name"] == "protac_3"][0]

        cham_flex = np.mean([r["flexibility_index"] for r in cham])
        cham_width = np.mean([r["spectral_width"] for r in cham])

        lines.append(f"\nChameleonic avg (protac_1/2):")
        lines.append(f"  Flexibility index: {cham_flex:.6f}")
        lines.append(f"  Spectral width: {cham_width:.4f}")

        lines.append(f"\nNon-chameleonic (protac_3):")
        lines.append(f"  Flexibility index: {non_cham['flexibility_index']:.6f}")
        lines.append(f"  Spectral width: {non_cham['spectral_width']:.4f}")

        diff = cham_flex - non_cham["flexibility_index"]
        lines.append(f"\nDifference: {diff:+.6f}")

        if diff < 0:
            lines.append(
                "\n>> SUCCESS: Non-chameleonic protac_3 shows HIGHER flexibility"
            )
        else:
            lines.append("\n>> INCONCLUSIVE: No clear separation")

    # Benchmark analysis
    lines.append("")
    lines.append("=" * 50)
    lines.append("BENCHMARK ANALYSIS")
    lines.append("=" * 50)

    benchmark_data = []
    descriptors = []

    try:
        with open("chameleon_local/benchmark.tsv", "r") as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    name, smiles, label = parts[0], parts[1], parts[2]
                    props = compute_graph_spectral_properties(smiles)
                    if props:
                        benchmark_data.append({"name": name, "label": label, **props})
                        descriptors.append(
                            {
                                "name": name,
                                "spectral_width": props["spectral_width"],
                                "spectral_gap": props["spectral_gap"],
                                "flexibility_index": props["flexibility_index"],
                                "label": label,
                            }
                        )
    except Exception as e:
        lines.append(f"Error loading benchmark: {e}")

    # Compute stats
    cham = [r for r in benchmark_data if r["label"] == "chameleonic"]
    non = [r for r in benchmark_data if r["label"] == "non_chameleonic"]

    if cham and non:
        lines.append(f"\nChameleonic (n={len(cham)}):")
        lines.append(
            f"  Flexibility index: {np.mean([r['flexibility_index'] for r in cham]):.6f}"
        )
        lines.append(
            f"  Spectral width: {np.mean([r['spectral_width'] for r in cham]):.4f}"
        )

        lines.append(f"\nNon-chameleonic (n={len(non)}):")
        lines.append(
            f"  Flexibility index: {np.mean([r['flexibility_index'] for r in non]):.6f}"
        )
        lines.append(
            f"  Spectral width: {np.mean([r['spectral_width'] for r in non]):.4f}"
        )

        diff = np.mean([r["flexibility_index"] for r in cham]) - np.mean(
            [r["flexibility_index"] for r in non]
        )
        lines.append(f"\nDifference: {diff:+.6f}")

    # Save descriptors
    with open("experiments/iter_24_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)
    lines.append("\nDescriptors saved to experiments/iter_24_descriptors.json")

    # Physical interpretation
    lines.append("")
    lines.append("=" * 50)
    lines.append("PHYSICAL INTERPRETATION")
    lines.append("=" * 50)
    lines.append("""
This method borrows from normal mode analysis in protein dynamics.
The spectral width measures the range of vibrational modes.

- High spectral width + high flexibility index: Very floppy, continuous motion
- Moderate values: Structured flexibility enabling discrete conformational states
- Low values: Rigid, limited flexibility

Chameleonic molecules should show MODERATE values - enough flexibility
to access distinct states, but not so floppy that they never stay folded.
""")

    # Final output
    output = "\n".join(lines)

    # Write to file
    with open("experiments/iter_24_output.txt", "w") as f:
        f.write(output)

    # Print to stderr
    print(output, file=sys.stderr)


if __name__ == "__main__":
    main()
