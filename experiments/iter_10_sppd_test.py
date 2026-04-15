
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))

def get_path(mol, start, end):
    # Simple BFS to find the shortest path
    q = [(start, [start])]
    visited = {start}
    while q:
        curr, path = q.pop(0)
        if curr == end:
            return path
        for neighbor in mol.GetAtomWithIdx(curr).GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited:
                visited.add(n_idx)
                q.append((n_idx, path + [n_idx]))
    return None

def get_shortest_path_polar_density(mol):
    try:
        dist_matrix = Chem.GetDistanceMatrix(mol)
        max_dist = np.max(dist_matrix)
        pairs = np.argwhere(dist_matrix == max_dist)
        
        densities = []
        for i, j in pairs:
            if i >= j: continue
            path = get_path(mol, int(i), int(j))
            if not path: continue
            
            polar_count = 0
            for idx in path:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() in ("N", "O"):
                    polar_count += 1
            densities.append(polar_count / len(path))
        
        if not densities: return 0.0
        return np.mean(densities)
    except Exception as e:
        print(f"Error in SPPD calculation: {e}")
        return 0.0

def run_exploration():
    try:
        user_mols = []
        with open("chameleon_local/user_protacs.tsv") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    user_mols.append((parts[0], parts[1]))

        print(f"{'Name':<20} | {'MaxDist':>7} | {'SPPD':>7}")
        print("-" * 40)

        for name, smi in user_mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(f"Failed to parse SMILES for {name}")
                continue
            max_dist = np.max(Chem.GetDistanceMatrix(mol))
            sppd = get_shortest_path_polar_density(mol)
            print(f"{name:<20} | {max_dist:>7.0f} | {sppd:>7.3f}")
    except Exception as e:
        print(f"Fatal error: {e}")

if __name__ == "__main__":
    run_exploration()
