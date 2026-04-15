import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

print("Starting iteration 27...")

from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import defaultdict
import math
import json


def load_molecules(filepath):
    print(f"Loading from {filepath}...")
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
    print(f"Loaded {len(molecules)} molecules")
    return molecules


def compute_path_entropy(mol):
    polar_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ["N", "O"]:
            polar_atoms.append(atom.GetIdx())

    if len(polar_atoms) < 2:
        return 0.0

    adj = defaultdict(list)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[i].append(j)
        adj[j].append(i)

    path_lengths = []
    for i, atom_i in enumerate(polar_atoms):
        for atom_j in polar_atoms[i + 1 :]:
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
        return 0.0

    path_counts = defaultdict(int)
    for pl in path_lengths:
        path_counts[pl] += 1

    total = len(path_lengths)
    entropy = 0.0

    for count in path_counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


def compute_linker_entropy(mol):
    terminal_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetDegree() == 1 and atom.GetSymbol() != "H":
            terminal_atoms.append(atom.GetIdx())

    if len(terminal_atoms) < 2:
        return 0.0

    adj = defaultdict(list)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj[i].append(j)
        adj[j].append(i)

    path_lengths = []
    for i, atom_i in enumerate(terminal_atoms):
        for atom_j in terminal_atoms[i + 1 :]:
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
        return 0.0

    path_counts = defaultdict(int)
    for pl in path_lengths:
        path_counts[pl] += 1

    total = len(path_lengths)
    entropy = 0.0

    for count in path_counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


def main():
    print("Loading molecules...")
    user_protacs = load_molecules("chameleon_local/user_protacs.tsv")
    benchmark = load_molecules("chameleon_local/benchmark.tsv")

    results = []
    output_lines = []

    output_lines.append("ITERATION 27: TOPOLOGICAL PATH ENTROPY")
    output_lines.append("=" * 60)
    output_lines.append("")
    output_lines.append("PROTAC RESULTS:")
    output_lines.append(f"{'Name':<12} {'PathEnt':<10} {'LinkEnt':<10}")
    output_lines.append("-" * 32)

    for name, smiles, mol in user_protacs:
        path_ent = compute_path_entropy(mol)
        link_ent = compute_linker_entropy(mol)

        results.append(
            {
                "name": name,
                "path_entropy": round(path_ent, 4),
                "linker_entropy": round(link_ent, 4),
            }
        )

        output_lines.append(f"{name:<12} {path_ent:<10.3f} {link_ent:<10.3f}")

    output_lines.append("")
    output_lines.append("Processing benchmark...")

    for name, smiles, mol in benchmark:
        path_ent = compute_path_entropy(mol)
        link_ent = compute_linker_entropy(mol)

        results.append(
            {
                "name": name,
                "path_entropy": round(path_ent, 4),
                "linker_entropy": round(link_ent, 4),
            }
        )

    output_lines.append(f"Processed {len(benchmark)} benchmark molecules")
    output_lines.append("")

    # Summary
    protac_results = [r for r in results if r["name"].startswith("protac_")]
    if len(protac_results) >= 3:
        protac_3 = [r for r in protac_results if r["name"] == "protac_3"][0]
        protac_12 = [r for r in protac_results if r["name"] in ["protac_1", "protac_2"]]

        if protac_12:
            avg_12_path = sum(r["path_entropy"] for r in protac_12) / len(protac_12)
            avg_12_link = sum(r["linker_entropy"] for r in protac_12) / len(protac_12)

            output_lines.append("SUMMARY:")
            output_lines.append(f"protac_1/2 avg path entropy: {avg_12_path:.3f}")
            output_lines.append(
                f"protac_3 path entropy: {protac_3['path_entropy']:.3f}"
            )
            output_lines.append(
                f"Path entropy gap: {avg_12_path - protac_3['path_entropy']:.3f}"
            )
            output_lines.append("")
            output_lines.append(f"protac_1/2 avg linker entropy: {avg_12_link:.3f}")
            output_lines.append(
                f"protac_3 linker entropy: {protac_3['linker_entropy']:.3f}"
            )
            output_lines.append(
                f"Linker entropy gap: {avg_12_link - protac_3['linker_entropy']:.3f}"
            )

    output_text = "\n".join(output_lines)

    with open("experiments/iter_27_output.txt", "w") as f:
        f.write(output_text)

    with open("experiments/iter_27_descriptors.json", "w") as f:
        json.dump(results, f, indent=2)

    print(output_text)


if __name__ == "__main__":
    main()
