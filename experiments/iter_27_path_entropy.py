import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import json
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import Descriptors
import math


def load_molecules(filepath):
    """Load molecules from TSV file."""
    molecules = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                name = parts[0]
                smiles = parts[1]
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    molecules.append((name, smiles, mol))
    return molecules


def compute_path_entropy(mol):
    """
    Compute Shannon entropy of shortest path lengths between polar atoms (N, O).
    Higher entropy = more diverse paths = chameleonic.

    Returns: (path_entropy, num_polar_pairs, unique_paths)
    """
    # Get all polar atoms (N, O)
    polar_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ["N", "O"]:
            polar_atoms.append(atom.GetIdx())

    if len(polar_atoms) < 2:
        return 0.0, 0, 0

    # Build adjacency list
    adj = defaultdict(list)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[i].append(j)
        adj[j].append(i)

    # Compute shortest paths between all polar atom pairs using BFS
    path_lengths = []

    for i, atom_i in enumerate(polar_atoms):
        for atom_j in polar_atoms[i + 1 :]:
            # BFS
            visited = {atom_i: 0}
            queue = [atom_i]
            found = False

            while queue and not found:
                current = queue.pop(0)
                current_dist = visited[current]

                for neighbor in adj[current]:
                    if neighbor not in visited:
                        visited[neighbor] = current_dist + 1
                        if neighbor == atom_j:
                            found = True
                            break
                        queue.append(neighbor)

            if atom_j in visited:
                path_lengths.append(visited[atom_j])

    if not path_lengths:
        return 0.0, 0, 0

    # Compute entropy of path length distribution
    # Bin by path length and compute probability
    path_counts = defaultdict(int)
    for pl in path_lengths:
        path_counts[pl] += 1

    total = len(path_lengths)
    entropy = 0.0

    for count in path_counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)

    return entropy, len(path_lengths), len(path_counts)


def compute_linker_path_entropy(mol):
    """
    Compute path entropy specifically for linker region.
    Focus on paths between terminal atoms to capture linker diversity.

    Returns: (linker_entropy, terminal_pairs)
    """
    # Find terminal atoms (degree 1, not hydrogen)
    terminal_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetDegree() == 1 and atom.GetSymbol() != "H":
            terminal_atoms.append(atom.GetIdx())

    if len(terminal_atoms) < 2:
        return 0.0, 0

    # Build adjacency
    adj = defaultdict(list)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[i].append(j)
        adj[j].append(i)

    # Compute paths between terminals
    path_lengths = []

    for i, atom_i in enumerate(terminal_atoms):
        for atom_j in terminal_atoms[i + 1 :]:
            # BFS
            visited = {atom_i: 0}
            queue = [atom_i]
            found = False

            while queue and not found:
                current = queue.pop(0)
                current_dist = visited[current]

                for neighbor in adj[current]:
                    if neighbor not in visited:
                        visited[neighbor] = current_dist + 1
                        if neighbor == atom_j:
                            found = True
                            break
                        queue.append(neighbor)

            if atom_j in visited:
                path_lengths.append(visited[atom_j])

    if not path_lengths:
        return 0.0, len(terminal_atoms)

    # Compute entropy
    path_counts = defaultdict(int)
    for pl in path_lengths:
        path_counts[pl] += 1

    total = len(path_lengths)
    entropy = 0.0

    for count in path_counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)

    return entropy, len(terminal_atoms)


def compute_weiner_index(mol):
    """
    Wiener index: sum of shortest path lengths between all atom pairs.
    Lower = more compact/connected graph.
    Higher = more extended/linear graph.
    """
    n_atoms = mol.GetNumAtoms()

    # Build adjacency
    adj = defaultdict(list)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[i].append(j)
        adj[j].append(i)

    # Compute all-pairs shortest paths
    total_distance = 0
    num_pairs = 0

    for start in range(n_atoms):
        # BFS from start
        visited = {start: 0}
        queue = [start]

        while queue:
            current = queue.pop(0)
            current_dist = visited[current]

            for neighbor in adj[current]:
                if neighbor not in visited:
                    visited[neighbor] = current_dist + 1
                    queue.append(neighbor)

        # Sum distances to all other atoms
        for end in range(start + 1, n_atoms):
            if end in visited:
                total_distance += visited[end]
                num_pairs += 1

    return total_distance, num_pairs


def compute_graph_spectral_width(mol):
    """
    Compute spectral width from graph Laplacian eigenvalues.
    Spectral width = ratio of largest to smallest non-zero eigenvalue.
    Measures graph connectivity and "stiffness".

    Returns: (spectral_width, num_eigenvalues)
    """
    import numpy as np

    n = mol.GetNumAtoms()
    if n < 3:
        return 0.0, 0

    # Build adjacency matrix
    adj = np.zeros((n, n))
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[i, j] = 1
        adj[j, i] = 1

    # Degree matrix
    degrees = np.sum(adj, axis=1)
    D = np.diag(degrees)

    # Laplacian
    L = D - adj

    # Eigenvalues
    try:
        eigenvalues = np.linalg.eigvalsh(L)
        eigenvalues = np.sort(eigenvalues)

        # Non-zero eigenvalues
        non_zero = eigenvalues[eigenvalues > 1e-10]

        if len(non_zero) < 2:
            return 0.0, len(eigenvalues)

        spectral_width = non_zero[-1] / non_zero[0]  # Max / Min
        return spectral_width, len(eigenvalues)
    except:
        return 0.0, 0


def main():
    # Load data
    user_protacs = load_molecules("chameleon_local/user_protacs.tsv")
    benchmark = load_molecules("chameleon_local/benchmark.tsv")

    all_molecules = user_protacs + benchmark

    results = []
    output_lines = []

    output_lines.append("=" * 70)
    output_lines.append("ITERATION 27: TOPOLOGICAL PATH ENTROPY")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append(
        "Hypothesis: PEG linkers create diverse path entropy (chameleonic)"
    )
    output_lines.append(
        "            Alkyl linkers create linear paths (non-chameleonic)"
    )
    output_lines.append("")
    output_lines.append("-" * 70)

    # Process PROTACs first
    output_lines.append("\nPROTAC RESULTS:")
    output_lines.append(
        f"{'Name':<12} {'Path_Ent':<10} {'Lnk_Ent':<10} {'Wiener':<10} {'SpectW':<10}"
    )
    output_lines.append("-" * 52)

    for name, smiles, mol in user_protacs:
        path_ent, n_pairs, n_unique = compute_path_entropy(mol)
        link_ent, n_term = compute_linker_path_entropy(mol)
        wiener, n_p = compute_weiner_index(mol)
        spect_w, n_eig = compute_graph_spectral_width(mol)
        mw = Descriptors.MolWt(mol)

        result = {
            "name": name,
            "path_entropy": round(path_ent, 4),
            "linker_entropy": round(link_ent, 4),
            "wiener_index": wiener,
            "spectral_width": round(spect_w, 4),
            "mw": mw,
        }
        results.append(result)

        output_lines.append(
            f"{name:<12} {path_ent:<10.3f} {link_ent:<10.3f} {wiener:<10} {spect_w:<10.3f}"
        )

    output_lines.append("")
    output_lines.append("-" * 70)
    output_lines.append("\nBENCHMARK MOLECULES (first 10):")
    output_lines.append(
        f"{'Name':<15} {'Path_Ent':<10} {'Lnk_Ent':<10} {'Wiener':<10} {'SpectW':<10}"
    )
    output_lines.append("-" * 55)

    for name, smiles, mol in benchmark[:10]:
        path_ent, n_pairs, n_unique = compute_path_entropy(mol)
        link_ent, n_term = compute_linker_path_entropy(mol)
        wiener, n_p = compute_weiner_index(mol)
        spect_w, n_eig = compute_graph_spectral_width(mol)
        mw = Descriptors.MolWt(mol)

        result = {
            "name": name,
            "path_entropy": round(path_ent, 4),
            "linker_entropy": round(link_ent, 4),
            "wiener_index": wiener,
            "spectral_width": round(spect_w, 4),
            "mw": mw,
        }
        results.append(result)

        output_lines.append(
            f"{name:<15} {path_ent:<10.3f} {link_ent:<10.3f} {wiener:<10} {spect_w:<10.3f}"
        )

    # Complete benchmark
    output_lines.append("")
    output_lines.append("Computing remaining benchmark molecules...")

    for name, smiles, mol in benchmark[10:]:
        path_ent, n_pairs, n_unique = compute_path_entropy(mol)
        link_ent, n_term = compute_linker_path_entropy(mol)
        wiener, n_p = compute_weiner_index(mol)
        spect_w, n_eig = compute_graph_spectral_width(mol)
        mw = Descriptors.MolWt(mol)

        result = {
            "name": name,
            "path_entropy": round(path_ent, 4),
            "linker_entropy": round(link_ent, 4),
            "wiener_index": wiener,
            "spectral_width": round(spect_w, 4),
            "mw": mw,
        }
        results.append(result)

    # Summary stats
    output_lines.append("")
    output_lines.append("=" * 70)
    output_lines.append("SUMMARY STATISTICS")
    output_lines.append("=" * 70)

    # PROTAC stats
    protac_results = [r for r in results if r["name"].startswith("protac_")]
    if len(protac_results) >= 3:
        path_ents = [r["path_entropy"] for r in protac_results]
        link_ents = [r["linker_entropy"] for r in protac_results]

        output_lines.append(
            f"\nPath Entropy Range: {min(path_ents):.3f} - {max(path_ents):.3f}"
        )
        output_lines.append(
            f"Linker Entropy Range: {min(link_ents):.3f} - {max(link_ents):.3f}"
        )

        # Check separation
        protac_3 = [r for r in protac_results if r["name"] == "protac_3"][0]
        protac_12 = [r for r in protac_results if r["name"] in ["protac_1", "protac_2"]]

        if protac_12:
            avg_12_path = sum(r["path_entropy"] for r in protac_12) / len(protac_12)
            avg_12_link = sum(r["linker_entropy"] for r in protac_12) / len(protac_12)

            output_lines.append(f"\nprotac_1/2 avg path entropy: {avg_12_path:.3f}")
            output_lines.append(
                f"protac_3 path entropy: {protac_3['path_entropy']:.3f}"
            )
            output_lines.append(
                f"Separation: {avg_12_path - protac_3['path_entropy']:.3f}"
            )

            output_lines.append(f"\nprotac_1/2 avg linker entropy: {avg_12_link:.3f}")
            output_lines.append(
                f"protac_3 linker entropy: {protac_3['linker_entropy']:.3f}"
            )
            output_lines.append(
                f"Separation: {avg_12_link - protac_3['linker_entropy']:.3f}"
            )

    output_lines.append("")
    output_lines.append("=" * 70)
    output_lines.append("PHYSICAL INTERPRETATION")
    output_lines.append("=" * 70)
    output_lines.append(
        "Path entropy measures the diversity of shortest paths between polar atoms."
    )
    output_lines.append(
        "PEG linkers (with ether oxygens) create multiple equivalent paths."
    )
    output_lines.append("Alkyl linkers create a single dominant path.")
    output_lines.append(
        "Higher entropy = more folding options = chameleonic propensity."
    )

    output_text = "\n".join(output_lines)

    # Write to file
    with open("experiments/iter_27_output.txt", "w") as f:
        f.write(output_text)

    # Also write descriptors
    with open("experiments/iter_27_descriptors.json", "w") as f:
        json.dump(results, f, indent=2)

    print(output_text)


if __name__ == "__main__":
    main()
