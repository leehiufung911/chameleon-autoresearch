#!/usr/bin/env python3
"""
Iteration 24: Molecular Spectral Flexibility Index
Graph spectral analysis from solid-state physics
"""

import os
import sys

os.environ["PYTHONUNBUFFERED"] = "1"

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import json

# Output collection
output_lines = []


def log(msg):
    """Log to both stdout and collection."""
    output_lines.append(str(msg))
    print(str(msg), flush=True)


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

    # Additional descriptors - simplified fallback
    n_rot = 0
    mw = 0.0
    try:
        from rdkit.Chem import Descriptors

        if hasattr(Descriptors, "NumRotatableBonds") and callable(
            getattr(Descriptors, "NumRotatableBonds", None)
        ):
            n_rot = getattr(Descriptors, "NumRotatableBonds")(mol)
        if hasattr(Descriptors, "MolWt") and callable(
            getattr(Descriptors, "MolWt", None)
        ):
            mw = getattr(Descriptors, "MolWt")(mol)
    except Exception:
        # Fallback: count rotatable bonds manually
        n_rot = 0
        for bond in mol.GetBonds():
            try:
                if bond.GetIsRotatable():
                    n_rot += 1
            except:
                pass
        mw = 0.0

    # Normalized metrics
    flex_index = spectral_width / (n_rot + 1)

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
    log("=" * 70)
    log("ITERATION 24: Molecular Spectral Flexibility Index")
    log("Graph spectral analysis from solid-state physics")
    log("=" * 70)
    log("")
    log(f"Working directory: {os.getcwd()}")
    log(f"Script location: {os.path.dirname(os.path.abspath(__file__))}")

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
    log("USER PROTAC ANALYSIS")
    log("-" * 50)

    user_results = []
    for name, smiles in user_protacs:
        props = compute_graph_spectral_properties(smiles)
        if props:
            user_results.append({"name": name, **props})
            log(f"\n{name}:")
            log(f"  Spectral width: {props['spectral_width']:.4f}")
            log(f"  Spectral gap: {props['spectral_gap']:.4f}")
            log(f"  Flexibility index: {props['flexibility_index']:.6f}")
            log(f"  Rotatable bonds: {props['n_rot_bonds']}")

    # Compare
    if len(user_results) == 3:
        log("")
        log("=" * 50)
        log("PROTAC COMPARISON")
        log("=" * 50)

        cham = [r for r in user_results if r["name"] in ["protac_1", "protac_2"]]
        non_cham = [r for r in user_results if r["name"] == "protac_3"][0]

        cham_flex = np.mean([r["flexibility_index"] for r in cham])
        cham_width = np.mean([r["spectral_width"] for r in cham])

        log(f"\nChameleonic avg (protac_1/2):")
        log(f"  Flexibility index: {cham_flex:.6f}")
        log(f"  Spectral width: {cham_width:.4f}")

        log(f"\nNon-chameleonic (protac_3):")
        log(f"  Flexibility index: {non_cham['flexibility_index']:.6f}")
        log(f"  Spectral width: {non_cham['spectral_width']:.4f}")

        diff = cham_flex - non_cham["flexibility_index"]
        log(f"\nDifference: {diff:+.6f}")

        if diff < 0:
            log("\n>> SUCCESS: Non-chameleonic protac_3 shows HIGHER flexibility")
        else:
            log("\n>> INCONCLUSIVE: No clear separation")

    # Benchmark analysis
    log("")
    log("=" * 50)
    log("BENCHMARK ANALYSIS")
    log("=" * 50)

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
        log(f"Error loading benchmark: {e}")

    # Compute stats
    cham = [r for r in benchmark_data if r["label"] == "chameleonic"]
    non = [r for r in benchmark_data if r["label"] == "non_chameleonic"]

    if cham and non:
        log(f"\nChameleonic (n={len(cham)}):")
        cham_flex_mean = np.mean([r["flexibility_index"] for r in cham])
        cham_width_mean = np.mean([r["spectral_width"] for r in cham])
        log(f"  Flexibility index: {cham_flex_mean:.6f}")
        log(f"  Spectral width: {cham_width_mean:.4f}")

        log(f"\nNon-chameleonic (n={len(non)}):")
        non_flex_mean = np.mean([r["flexibility_index"] for r in non])
        non_width_mean = np.mean([r["spectral_width"] for r in non])
        log(f"  Flexibility index: {non_flex_mean:.6f}")
        log(f"  Spectral width: {non_width_mean:.4f}")

        diff = cham_flex_mean - non_flex_mean
        log(f"\nDifference: {diff:+.6f}")

    # Save descriptors
    with open("experiments/iter_24_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)
    log("\nDescriptors saved to experiments/iter_24_descriptors.json")

    # Physical interpretation
    log("")
    log("=" * 50)
    log("PHYSICAL INTERPRETATION")
    log("=" * 50)
    log("This method borrows from normal mode analysis in protein dynamics.")
    log("The spectral width measures the range of vibrational modes.")
    log("")
    log("- High spectral width + high flexibility index: Very floppy")
    log("- Moderate values: Structured flexibility (chameleonic)")
    log("- Low values: Rigid, limited flexibility")
    log("")
    log("Chameleonic molecules should show MODERATE values - enough flexibility")
    log("to access distinct states, but not so floppy that they never stay folded.")

    # Write output file
    with open("experiments/iter_24_output.txt", "w") as f:
        f.write("\n".join(output_lines))

    log("\nOutput saved to experiments/iter_24_output.txt")


if __name__ == "__main__":
    main()
