
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))
import chameleon

def get_max_nonring_carbon_path(mol):
    nonring_c = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C" and not a.IsInRing()]
    if not nonring_c: return 0
    n = len(nonring_c)
    adj = np.zeros((n, n))
    idx_map = {idx: i for i, idx in enumerate(nonring_c)}
    for i, idx in enumerate(nonring_c):
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in idx_map:
                j = idx_map[neighbor.GetIdx()]
                adj[i, j] = 1
                adj[j, i] = 1
    dist = adj.copy()
    dist[dist == 0] = 999
    np.fill_diagonal(dist, 0)
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if dist[i, j] > dist[i, k] + dist[k, j]:
                    dist[i, j] = dist[i, k] + dist[k, j]
    dist[dist >= 999] = 0
    return int(np.max(dist)) + 1 if n > 0 else 0

def run_experiment():
    # Test on a few non-chameleons that are often false positives
    targets = ["atorvastatin", "warfarin", "metoprolol", "atenolol", "protac_3"]
    all_items = []
    with open("chameleon_local/benchmark.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] in targets:
                all_items.append((parts[0], parts[1]))
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] in targets:
                all_items.append((parts[0], parts[1]))

    print(f"{'Name':<20} | {'MaxC':>4} | {'IMHB':>5}")
    print("-" * 35)

    for name, smi in all_items:
        try:
            mol = Chem.MolFromSmiles(smi)
            max_c = get_max_nonring_carbon_path(mol)
            summary = chameleon.summarize(name, smi, n_conf=20, verbose=False)
            print(f"{name:<20} | {max_c:>4} | {summary.imhb_mean:>5.2f}")
        except Exception as e:
            pass

if __name__ == "__main__":
    run_experiment()
