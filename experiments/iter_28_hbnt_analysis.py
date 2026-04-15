import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from collections import defaultdict
import json


def get_hbond_donors_acceptor_indices(mol):
    """Get indices of H-bond donors and acceptors."""
    donors = []
    acceptors = []

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        atomic_num = atom.GetAtomicNum()

        # H-bond donors: N-H, O-H
        if atomic_num in [7, 8]:  # N or O
            # Check if it has an H attached
            has_h = any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors())
            if has_h:
                donors.append(idx)
            # If it has lone pairs, it's also an acceptor
            if len(atom.GetNeighbors()) < 4:
                acceptors.append(idx)

    return donors, acceptors


def compute_hbond_network_topology(smiles, name, n_conformers=20):
    """
    Analyze hydrogen bond network topology across conformer ensemble.
    Returns metrics on cooperative H-bond networks.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    # Generate conformers
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0

    AllChem.EmbedMultipleConfs(mol, n_conformers, params)

    if mol.GetNumConformers() == 0:
        return None

    # Get H-bond donors/acceptors
    donors, acceptors = get_hbond_donors_acceptor_indices(mol)

    if len(donors) == 0 or len(acceptors) == 0:
        return None

    hbond_stats = []

    for conf_id in range(mol.GetNumConformers()):
        conf = mol.GetConformer(conf_id)

        # Find all IMHBs in this conformer
        hbonds = []
        for d_idx in donors:
            d_atom = mol.GetAtomWithIdx(d_idx)
            d_pos = conf.GetAtomPosition(d_idx)

            for a_idx in acceptors:
                if d_idx == a_idx:
                    continue

                a_atom = mol.GetAtomWithIdx(a_idx)
                a_pos = conf.GetAtomPosition(a_idx)

                # Compute distance
                dist = np.sqrt(
                    (d_pos.x - a_pos.x) ** 2
                    + (d_pos.y - a_pos.y) ** 2
                    + (d_pos.z - a_pos.z) ** 2
                )

                # H-bond if distance < 2.5 Angstrom
                if dist < 2.5:
                    hbonds.append((d_idx, a_idx, dist))

        # Analyze network topology
        if len(hbonds) > 0:
            # Build adjacency list
            hbond_graph = defaultdict(list)
            for d, a, dist in hbonds:
                hbond_graph[d].append((a, dist))
                hbond_graph[a].append((d, dist))

            # Compute connected components (cooperative networks)
            visited = set()
            components = []

            for node in hbond_graph:
                if node not in visited:
                    component = []
                    stack = [node]
                    while stack:
                        current = stack.pop()
                        if current not in visited:
                            visited.add(current)
                            component.append(current)
                            for neighbor, _ in hbond_graph[current]:
                                if neighbor not in visited:
                                    stack.append(neighbor)
                    components.append(component)

            # Calculate network metrics
            max_component_size = max(len(c) for c in components) if components else 0
            num_components = len(components)
            avg_component_size = (
                np.mean([len(c) for c in components]) if components else 0
            )

            # Density of H-bond network
            network_density = (
                len(hbonds) / (len(donors) * len(acceptors))
                if (len(donors) * len(acceptors)) > 0
                else 0
            )

            # Check if there's a spanning network (connects across molecule)
            # Heuristic: largest component > 50% of H-bond capable atoms
            spanning = max_component_size > (len(donors) + len(acceptors)) * 0.5

            hbond_stats.append(
                {
                    "n_hbonds": len(hbonds),
                    "max_component_size": max_component_size,
                    "num_components": num_components,
                    "avg_component_size": avg_component_size,
                    "network_density": network_density,
                    "spanning": spanning,
                }
            )

    if not hbond_stats:
        return None

    # Aggregate across conformers
    mean_hbonds = np.mean([s["n_hbonds"] for s in hbond_stats])
    mean_max_component = np.mean([s["max_component_size"] for s in hbond_stats])
    mean_components = np.mean([s["num_components"] for s in hbond_stats])
    spanning_fraction = np.mean([1 if s["spanning"] else 0 for s in hbond_stats])
    mean_density = np.mean([s["network_density"] for s in hbond_stats])

    # Cooperativity metric: large components vs many small ones
    # High = cooperative, Low = fragmented
    cooperativity = mean_max_component / (mean_components + 1)

    return {
        "name": name,
        "n_donors": len(donors),
        "n_acceptors": len(acceptors),
        "mean_hbonds": mean_hbonds,
        "mean_max_component": mean_max_component,
        "mean_components": mean_components,
        "spanning_fraction": spanning_fraction,
        "network_density": mean_density,
        "cooperativity": cooperativity,
    }


def main():
    output = []
    output.append("=" * 70)
    output.append("ITERATION 28: Hydrogen Bond Network Topology Analysis")
    output.append("=" * 70)
    output.append("")
    output.append("Hypothesis: Chameleonic molecules form cooperative H-bond networks")
    output.append("that span across the molecule (high betweenness, large connected")
    output.append("components), while non-chameleonic molecules have fragmented,")
    output.append("non-cooperative H-bonds.")
    output.append("")
    output.append("Physical Mechanism:")
    output.append("- Chameleonic: Coordinated folding creates spanning H-bond networks")
    output.append("- Non-chameleonic: Local H-bonds only, no cooperative stabilization")
    output.append("- Metric: max_component_size / num_components = cooperativity")
    output.append("")

    # Load PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) >= 2:
                name = parts[0]
                smiles = parts[1].strip()
                protacs.append((name, smiles))

    labels = {"protac_1": "CHAM", "protac_2": "CHAM", "protac_3": "NON"}

    results = []
    output.append("Processing PROTACs...")
    output.append("-" * 70)

    for name, smiles in protacs:
        output.append(f"\nProcessing {name}...")
        result = compute_hbond_network_topology(smiles, name, n_conformers=15)

        if result:
            results.append(result)
            label = labels.get(name, "?")
            output.append(f"  Label: {label}")
            output.append(
                f"  H-bond donors: {result['n_donors']}, acceptors: {result['n_acceptors']}"
            )
            output.append(f"  Mean H-bonds: {result['mean_hbonds']:.2f}")
            output.append(f"  Max component size: {result['mean_max_component']:.2f}")
            output.append(f"  Num components: {result['mean_components']:.2f}")
            output.append(f"  Spanning fraction: {result['spanning_fraction']:.2f}")
            output.append(f"  Network density: {result['network_density']:.4f}")
            output.append(f"  **Cooperativity**: {result['cooperativity']:.3f}")
        else:
            output.append(f"  FAILED to analyze")

    # Summary
    output.append("\n" + "=" * 70)
    output.append("SUMMARY")
    output.append("=" * 70)

    if len(results) >= 3:
        cham_coop = [
            r["cooperativity"] for r in results if labels.get(r["name"]) == "CHAM"
        ]
        non_coop = [
            r["cooperativity"] for r in results if labels.get(r["name"]) == "NON"
        ]

        if cham_coop and non_coop:
            cham_mean = np.mean(cham_coop)
            non_mean = np.mean(non_coop)
            gap = cham_mean - non_mean

            output.append(f"\nChameleonic avg cooperativity: {cham_mean:.3f}")
            output.append(f"Non-chameleonic avg: {non_mean:.3f}")
            output.append(f"Separation: {gap:+.3f} (positive = correct)")

            output.append("\nResults table:")
            output.append(
                f"{'Name':<12} {'Label':<6} {'Coop':<8} {'MaxComp':<8} {'SpanFrac':<8}"
            )
            output.append("-" * 50)
            for r in results:
                label = labels.get(r["name"], "?")
                output.append(
                    f"{r['name']:<12} {label:<6} {r['cooperativity']:.3f}    "
                    f"{r['mean_max_component']:.2f}     {r['spanning_fraction']:.2f}"
                )

            if gap > 0.05:
                output.append(
                    "\nVerdict: SUCCESS - Chameleonic PROTACs show higher H-bond cooperativity"
                )
            elif gap > 0:
                output.append(
                    "\nVerdict: PARTIAL SUCCESS - Correct ordering but small gap"
                )
            else:
                output.append("\nVerdict: FAILED - Wrong ordering")

    # Write output
    output_text = "\n".join(output)
    with open("experiments/iter_28_output.txt", "w") as f:
        f.write(output_text)

    # Write descriptors JSON
    descriptors = []
    for r in results:
        descriptors.append(
            {
                "name": r["name"],
                "cooperativity": r["cooperativity"],
                "mean_max_component": r["mean_max_component"],
                "spanning_fraction": r["spanning_fraction"],
                "network_density": r["network_density"],
            }
        )

    with open("experiments/iter_28_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)

    print(output_text)
    print("\nOutput written to experiments/iter_28_output.txt")
    print("Descriptors written to experiments/iter_28_descriptors.json")


if __name__ == "__main__":
    main()
