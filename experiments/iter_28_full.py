import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import json


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
            from collections import defaultdict

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


# Load benchmark molecules
results = []
with open("chameleon_local/benchmark.tsv", "r") as f:
    for line in f:
        parts = line.strip().split(maxsplit=1)
        if len(parts) >= 2:
            name = parts[0]
            smiles = parts[1].strip()
            result = compute_hbond_network_topology(smiles, name, n_conformers=15)
            if result:
                results.append(result)
                sys.stderr.write(
                    f"{name}: cooperativity = {result['cooperativity']:.3f}\n"
                )
                sys.stderr.flush()

# Load PROTACs
protac_results = []
with open("chameleon_local/user_protacs.tsv", "r") as f:
    for line in f:
        parts = line.strip().split(maxsplit=1)
        if len(parts) >= 2:
            name = parts[0]
            smiles = parts[1].strip()
            result = compute_hbond_network_topology(smiles, name, n_conformers=15)
            if result:
                protac_results.append(result)
                sys.stderr.write(
                    f"{name}: cooperativity = {result['cooperativity']:.3f}\n"
                )
                sys.stderr.flush()

# Write descriptors for all molecules
all_results = results + protac_results
with open("experiments/iter_28_descriptors.json", "w") as f:
    json.dump(all_results, f, indent=2)

sys.stderr.write(
    f"\nWrote {len(all_results)} descriptors to iter_28_descriptors.json\n"
)
sys.stderr.flush()
