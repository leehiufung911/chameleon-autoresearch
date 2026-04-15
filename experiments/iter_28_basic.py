import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import json
from collections import defaultdict

print("Script starting")


def get_hbond_donors_acceptor_indices(mol):
    donors = []
    acceptors = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        atomic_num = atom.GetAtomicNum()
        if atomic_num in [7, 8]:
            has_h = any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors())
            if has_h:
                donors.append(idx)
            if len(atom.GetNeighbors()) < 4:
                acceptors.append(idx)
    return donors, acceptors


def compute_hbond_network_topology(smiles, name, n_conformers=15):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0

    AllChem.EmbedMultipleConfs(mol, n_conformers, params)
    if mol.GetNumConformers() == 0:
        return None

    donors, acceptors = get_hbond_donors_acceptor_indices(mol)
    if len(donors) == 0 or len(acceptors) == 0:
        return None

    hbond_stats = []

    for conf_id in range(mol.GetNumConformers()):
        conf = mol.GetConformer(conf_id)
        hbonds = []
        for d_idx in donors:
            d_pos = conf.GetAtomPosition(d_idx)
            for a_idx in acceptors:
                if d_idx == a_idx:
                    continue
                a_pos = conf.GetAtomPosition(a_idx)
                dist = np.sqrt(
                    (d_pos.x - a_pos.x) ** 2
                    + (d_pos.y - a_pos.y) ** 2
                    + (d_pos.z - a_pos.z) ** 2
                )
                if dist < 2.5:
                    hbonds.append((d_idx, a_idx, dist))

        if len(hbonds) > 0:
            hbond_graph = defaultdict(list)
            for d, a, dist in hbonds:
                hbond_graph[d].append((a, dist))
                hbond_graph[a].append((d, dist))

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

            max_component_size = max(len(c) for c in components) if components else 0
            num_components = len(components)
            spanning = max_component_size > (len(donors) + len(acceptors)) * 0.5

            hbond_stats.append(
                {
                    "n_hbonds": len(hbonds),
                    "max_component_size": max_component_size,
                    "num_components": num_components,
                    "spanning": spanning,
                }
            )

    if not hbond_stats:
        return None

    mean_max_component = np.mean([s["max_component_size"] for s in hbond_stats])
    mean_components = np.mean([s["num_components"] for s in hbond_stats])
    spanning_fraction = np.mean([1 if s["spanning"] else 0 for s in hbond_stats])
    cooperativity = mean_max_component / (mean_components + 1)

    return {
        "name": name,
        "mean_max_component": float(mean_max_component),
        "mean_components": float(mean_components),
        "spanning_fraction": float(spanning_fraction),
        "cooperativity": float(cooperativity),
    }


print("Functions defined")

# Load PROTACs
protacs = []
with open("chameleon_local/user_protacs.tsv", "r") as f:
    for line in f:
        parts = line.strip().split(maxsplit=1)
        if len(parts) >= 2:
            name = parts[0]
            smiles = parts[1].strip()
            protacs.append((name, smiles))

print(f"Loaded {len(protacs)} PROTACs")

labels = {"protac_1": "CHAM", "protac_2": "CHAM", "protac_3": "NON"}

results = []

for name, smiles in protacs:
    print(f"Processing {name}...")
    result = compute_hbond_network_topology(smiles, name, n_conformers=15)
    if result:
        results.append(result)
        label = labels.get(name, "?")
        print(f"  {label}: Cooperativity = {result['cooperativity']:.3f}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

if len(results) >= 3:
    cham_coop = [r["cooperativity"] for r in results if labels.get(r["name"]) == "CHAM"]
    non_coop = [r["cooperativity"] for r in results if labels.get(r["name"]) == "NON"]

    if cham_coop and non_coop:
        cham_mean = np.mean(cham_coop)
        non_mean = np.mean(non_coop)
        gap = cham_mean - non_mean

        print(f"Chameleonic avg cooperativity: {cham_mean:.3f}")
        print(f"Non-chameleonic avg: {non_mean:.3f}")
        print(f"Gap: {gap:+.3f}")

        if gap > 0.05:
            print("Verdict: SUCCESS")
        elif gap > 0:
            print("Verdict: PARTIAL SUCCESS")
        else:
            print("Verdict: FAILED")

# Write output
output_lines = ["ITERATION 28: Hydrogen Bond Network Topology Analysis", ""]
for r in results:
    output_lines.append(f"{r['name']}: Cooperativity = {r['cooperativity']:.3f}")
output_lines.append("")
if len(results) >= 3:
    output_lines.append(f"Gap: {gap:+.3f}")

output_text = "\n".join(output_lines)

with open("experiments/iter_28_output.txt", "w") as f:
    f.write(output_text)

with open("experiments/iter_28_descriptors.json", "w") as f:
    json.dump(results, f, indent=2)

print("\nFiles written")
