import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import json
from collections import defaultdict


def compute_graph_laplacian(mol):
    """Compute graph Laplacian matrix for the molecule."""
    n_atoms = mol.GetNumAtoms()
    adjacency = np.zeros((n_atoms, n_atoms))

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        # Weight by bond type: single=1, double=2, aromatic=1.5
        if bond.GetIsAromatic():
            weight = 1.5
        else:
            bond_type = bond.GetBondType()
            weight = {
                Chem.BondType.SINGLE: 1.0,
                Chem.BondType.DOUBLE: 2.0,
                Chem.BondType.TRIPLE: 3.0,
            }.get(bond_type, 1.0)
        adjacency[i, j] = weight
        adjacency[j, i] = weight

    # Degree matrix
    degree = np.diag(adjacency.sum(axis=1))

    # Laplacian: L = D - A
    laplacian = degree - adjacency

    return laplacian


def compute_compliance_matrix(mol):
    """Compute approximate compliance matrix from graph structure."""
    n_atoms = mol.GetNumAtoms()

    # Build distance matrix (topological distances)
    distances = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            # Use topological distance as proxy
            try:
                path = Chem.rdmolops.GetShortestPath(mol, i, j)
                dist = len(path) - 1  # Number of bonds
            except:
                dist = n_atoms  # Disconnected
            distances[i, j] = dist
            distances[j, i] = dist

    # Build spring constant matrix
    # Closer atoms (in graph) have stronger springs
    k_matrix = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if distances[i, j] > 0:
                # Inverse square of distance as spring constant
                k = 1.0 / (distances[i, j] ** 2)
                k_matrix[i, j] = k
                k_matrix[j, i] = k

    # Build approximate Hessian (simplified)
    hessian = np.diag(k_matrix.sum(axis=1)) - k_matrix

    return hessian, distances


def compute_flexibility_index(mol):
    """
    Compute flexibility index from spectral properties.
    Lower index = more rigid/structured flexibility
    Higher index = more floppy/continuous flexibility
    Returns: (spectral_width, spectral_gap, effective_flexibility, hessian_width)
    """
    n_atoms = mol.GetNumAtoms()

    # Compute Laplacian
    laplacian = compute_graph_laplacian(mol)

    # Compute eigenvalues (spectral analysis)
    try:
        eigenvalues = np.linalg.eigvalsh(laplacian)
        # Remove zero eigenvalue (connected component)
        non_zero = eigenvalues[eigenvalues > 1e-10]

        if len(non_zero) == 0:
            return 0.0, 0.0, 0.0, 0.0

        # Spectral metrics
        smallest = np.min(non_zero)
        largest = np.max(non_zero)
        spectral_width = largest / smallest if smallest > 0 else 0

        # Spectral gap (difference between first and second smallest)
        sorted_vals = np.sort(non_zero)
        spectral_gap = sorted_vals[1] - sorted_vals[0] if len(sorted_vals) > 1 else 0

        # Effective stiffness (inverse of compliance)
        # Chameleonic molecules should have intermediate values
        # Not too stiff (small width), not too floppy (large width)
        effective_flexibility = spectral_width / n_atoms

    except Exception as e:
        print(f"Error computing eigenvalues: {e}")
        return 0.0, 0.0, 0.0, 0.0

    # Also compute graph-theoretic metrics
    hessian, distances = compute_compliance_matrix(mol)
    try:
        hessian_eig = np.linalg.eigvalsh(hessian)
        hessian_nonzero = hessian_eig[hessian_eig > 1e-10]
        if len(hessian_nonzero) > 0:
            hessian_width = np.max(hessian_nonzero) / np.min(hessian_nonzero)
        else:
            hessian_width = 0.0
    except:
        hessian_width = 0.0

    return spectral_width, spectral_gap, effective_flexibility, hessian_width


def analyze_molecule(name, smiles):
    """Analyze a single molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    mw = Descriptors.MolWt(mol)
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()

    # Compute flexibility metrics
    result = compute_flexibility_index(mol)
    spec_width = result[0]
    spec_gap = result[1]
    eff_flex = result[2]
    hess_width = result[3]

    # Compute additional descriptors
    rot_bonds = Descriptors.NumRotatableBonds(mol)

    return {
        "name": name,
        "mw": mw,
        "n_atoms": n_atoms,
        "n_bonds": n_bonds,
        "rotatable_bonds": rot_bonds,
        "spectral_width": spec_width,
        "spectral_gap": spec_gap,
        "effective_flexibility": eff_flex,
        "hessian_width": hess_width,
        "flexibility_per_atom": eff_flex / n_atoms if n_atoms > 0 else 0,
        "flexibility_per_rot": eff_flex / rot_bonds if rot_bonds > 0 else eff_flex,
    }


def main():
    output_lines = []

    # Load user PROTACs
    user_protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                name, smiles = parts
                user_protacs.append((name, smiles))

    output_lines.append("=" * 70)
    output_lines.append("ITERATION 24: Molecular Compliance (Flexibility Index)")
    output_lines.append("Graph Spectral Analysis from Solid-State Physics")
    output_lines.append("=" * 70)
    output_lines.append("")

    # Analyze user PROTACs
    output_lines.append("USER PROTAC ANALYSIS")
    output_lines.append("-" * 50)

    user_results = []
    for name, smiles in user_protacs:
        result = analyze_molecule(name, smiles)
        if result:
            user_results.append(result)
            output_lines.append(f"\n{name}:")
            output_lines.append(
                f"  MW: {result['mw']:.1f}, Atoms: {result['n_atoms']}, Bonds: {result['n_bonds']}"
            )
            output_lines.append(f"  Rotatable bonds: {result['rotatable_bonds']}")
            output_lines.append(f"  Spectral width: {result['spectral_width']:.4f}")
            output_lines.append(f"  Spectral gap: {result['spectral_gap']:.4f}")
            output_lines.append(
                f"  Effective flexibility: {result['effective_flexibility']:.6f}"
            )
            output_lines.append(
                f"  Flexibility per rotatable bond: {result['flexibility_per_rot']:.6f}"
            )
            output_lines.append(f"  Hessian width: {result['hessian_width']:.4f}")

    output_lines.append("")
    output_lines.append("=" * 50)
    output_lines.append("PROTAC COMPARISON")
    output_lines.append("=" * 50)

    # Compare key metrics
    chams = [r for r in user_results if r["name"] in ["protac_1", "protac_2"]]
    non_cham = [r for r in user_results if r["name"] == "protac_3"]

    if chams and non_cham:
        cham_spec = np.mean([r["spectral_width"] for r in chams])
        cham_flex = np.mean([r["flexibility_per_rot"] for r in chams])
        non_spec = non_cham[0]["spectral_width"]
        non_flex = non_cham[0]["flexibility_per_rot"]

        output_lines.append(f"\nChameleonic avg (protac_1/2):")
        output_lines.append(f"  Spectral width: {cham_spec:.4f}")
        output_lines.append(f"  Flexibility per rot: {cham_flex:.6f}")

        output_lines.append(f"\nNon-chameleonic (protac_3):")
        output_lines.append(f"  Spectral width: {non_spec:.4f}")
        output_lines.append(f"  Flexibility per rot: {non_flex:.6f}")

        output_lines.append(f"\nSeparation:")
        output_lines.append(f"  Spectral width: {cham_spec - non_spec:+.4f}")
        output_lines.append(f"  Flexibility per rot: {cham_flex - non_flex:+.6f}")

    # Load and analyze benchmark molecules
    output_lines.append("")
    output_lines.append("=" * 50)
    output_lines.append("BENCHMARK ANALYSIS (49 molecules)")
    output_lines.append("=" * 50)

    benchmark_results = []
    with open("chameleon_local/benchmark.tsv", "r") as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                name = parts[0]
                smiles = parts[1]
                label = parts[2]
                result = analyze_molecule(name, smiles)
                if result:
                    result["label"] = label
                    benchmark_results.append(result)

    # Compute statistics
    cham_results = [r for r in benchmark_results if r.get("label") == "chameleonic"]
    non_results = [r for r in benchmark_results if r.get("label") == "non_chameleonic"]

    if cham_results and non_results:
        cham_spec_mean = np.mean([r["spectral_width"] for r in cham_results])
        cham_flex_mean = np.mean([r["flexibility_per_rot"] for r in cham_results])
        non_spec_mean = np.mean([r["spectral_width"] for r in non_results])
        non_flex_mean = np.mean([r["flexibility_per_rot"] for r in non_results])

        output_lines.append(f"\nBenchmark Statistics:")
        output_lines.append(f"  Chameleonic (n={len(cham_results)}):")
        output_lines.append(f"    Avg spectral width: {cham_spec_mean:.4f}")
        output_lines.append(f"    Avg flexibility/rot: {cham_flex_mean:.6f}")
        output_lines.append(f"  Non-chameleonic (n={len(non_results)}):")
        output_lines.append(f"    Avg spectral width: {non_spec_mean:.4f}")
        output_lines.append(f"    Avg flexibility/rot: {non_flex_mean:.6f}")
        output_lines.append(f"  Separation:")
        output_lines.append(
            f"    Spectral width: {cham_spec_mean - non_spec_mean:+.4f}"
        )
        output_lines.append(
            f"    Flexibility/rot: {cham_flex_mean - non_flex_mean:+.6f}"
        )

    # Write descriptors to JSON
    all_descriptors = []
    for r in user_results + benchmark_results:
        desc = {
            "name": r["name"],
            "spectral_width": r["spectral_width"],
            "spectral_gap": r["spectral_gap"],
            "effective_flexibility": r["effective_flexibility"],
            "flexibility_per_rot": r["flexibility_per_rot"],
            "hessian_width": r["hessian_width"],
        }
        if "label" in r:
            desc["label"] = r["label"]
        all_descriptors.append(desc)

    with open("experiments/iter_24_descriptors.json", "w") as f:
        json.dump(all_descriptors, f, indent=2)

    output_lines.append("\nDescriptors written to experiments/iter_24_descriptors.json")

    output_lines.append("\n" + "=" * 50)
    output_lines.append("PHYSICAL INTERPRETATION")
    output_lines.append("=" * 50)
    output_lines.append("""
Spectral width measures the range of vibrational frequencies in the molecule.
- Low spectral width: Rigid structure, limited conformational freedom
- High spectral width: Very floppy, continuous flexibility
- Intermediate (chameleonic): Structured flexibility enabling discrete states

Flexibility per rotatable bond normalizes by the number of flexible joints.
- Chameleonic molecules should show LOWER values (structured flexibility)
- Non-chameleonic floppy molecules show HIGHER values (random flexibility)

The hypothesis predicts protac_3 will have HIGHER flexibility per rot bond
than protac_1/2, indicating continuous rather than discrete flexibility.
""")

    output_lines.append("\n" + "=" * 50)
    output_lines.append("CONCLUSION")
    output_lines.append("=" * 50)

    # Final assessment
    if cham_flex_mean and non_flex_mean:
        if non_flex_mean > cham_flex_mean:
            output_lines.append(
                "\n✓ SUCCESS: Non-chameleonic molecules show higher normalized flexibility"
            )
            output_lines.append(
                "  (continuous floppy motion vs structured discrete states)"
            )
        else:
            output_lines.append(
                "\n✗ INCONCLUSIVE: No clear separation in flexibility patterns"
            )

    output_text = "\n".join(output_lines)

    # Write to file
    with open("experiments/iter_24_output.txt", "w") as f:
        f.write(output_text)

    print(output_text)


if __name__ == "__main__":
    main()
