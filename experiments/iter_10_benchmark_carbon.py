
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))

def get_max_nonring_carbon_path(mol):
    nonring_c = []
    for a in mol.GetAtoms():
        if a.GetSymbol() == "C" and not a.IsInRing():
            nonring_c.append(a.GetIdx())
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
    mols = []
    with open("chameleon_local/labelled_set.tsv") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                mols.append((parts[0], parts[1], parts[2]))

    print(f"{'Label':<12} | {'Name':<20} | {'MaxC'}")
    print("-" * 45)

    for label, name, smi in mols:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            max_c = get_max_nonring_carbon_path(mol)
            if max_c >= 4:
                print(f"{label:<12} | {name:<20} | {max_c}")

if __name__ == "__main__":
    run_experiment()
