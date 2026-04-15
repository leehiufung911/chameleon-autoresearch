
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))
import chameleon

def get_linker_atoms(mol):
    # Find all rings
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return []
    
    # Group atoms by ring system
    ring_systems = []
    for ring in rings:
        ring_set = set(ring)
        # Merge with existing systems if overlapping
        merged = False
        for i, sys in enumerate(ring_systems):
            if ring_set.intersection(sys):
                ring_systems[i] = sys.union(ring_set)
                merged = True
                break
        if not merged:
            ring_systems.append(ring_set)
    
    if len(ring_systems) < 2:
        return []
        
    # Pick the two largest ring systems
    ring_systems.sort(key=len, reverse=True)
    sys1 = list(ring_systems[0])
    sys2 = list(ring_systems[1])
    
    # Find the shortest path between any atom in sys1 and any atom in sys2
    dist_matrix = Chem.GetDistanceMatrix(mol)
    min_dist = 999
    best_pair = (None, None)
    for a1 in sys1:
        for a2 in sys2:
            d = dist_matrix[a1, a2]
            if d < min_dist:
                min_dist = d
                best_pair = (a1, a2)
    
    if best_pair[0] is None:
        return []
        
    path = get_path(mol, int(best_pair[0]), int(best_pair[1]))
    # The "linker" atoms are the ones between the rings
    linker_atoms = [idx for idx in path if idx not in sys1 and idx not in sys2]
    return linker_atoms

def get_path(mol, start, end):
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

def compute_linker_imhb(mol, confId, linker_atoms, topo_dist):
    if not linker_atoms:
        return 0
    linker_set = set(linker_atoms)
    
    # Use chameleon's logic but filter for linker involvement
    # Actually, chameleon.compute_imhb returns the count.
    # I'll re-implement it briefly with linker filter.
    conf = mol.GetConformer(confId)
    count = 0
    
    # Find H-bond donors and acceptors
    donors = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in ("N", "O"):
            for n in a.GetNeighbors():
                if n.GetSymbol() == "H":
                    donors.append((a.GetIdx(), n.GetIdx()))
                    break
    acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ("N", "O")]
    
    for d_idx, h_idx in donors:
        for a_idx in acceptors:
            if d_idx == a_idx: continue
            if topo_dist[d_idx, a_idx] < 4: continue
            
            # Check if at least one atom is a linker atom
            if d_idx not in linker_set and a_idx not in linker_set:
                continue
                
            dist = np.linalg.norm(conf.GetAtomPosition(h_idx) - conf.GetAtomPosition(a_idx))
            if dist < 2.5: # Standard threshold
                count += 1
    return count

def run_experiment():
    user_mols = []
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                user_mols.append((parts[0], parts[1]))

    N_CONF = 20
    print(f"{'Name':<20} | {'IMHB_all':>8} | {'IMHB_link':>9} | {'% Linker'}")
    print("-" * 55)

    for name, smi in user_mols:
        try:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            linker_atoms = get_linker_atoms(mol)
            print(f"Linker atoms for {name}: {linker_atoms}")
            for idx in linker_atoms:
                atom = mol.GetAtomWithIdx(idx)
                print(f"  {atom.GetSymbol()}({idx})", end="")
            print()
            topo_dist = chameleon.bond_path_distance_matrix(mol)
            
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xC0FFEE
            cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=N_CONF, params=params))
            AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500, mmffVariant="MMFF94s")
            
            imhb_all = []
            imhb_link = []
            for cid in cids:
                imhb_all.append(chameleon.compute_imhb(mol, cid, topo_dist))
                imhb_link.append(compute_linker_imhb(mol, cid, linker_atoms, topo_dist))
            
            m_all = np.mean(imhb_all)
            m_link = np.mean(imhb_link)
            pct = (m_link / m_all * 100) if m_all > 0 else 0
            
            print(f"{name:<20} | {m_all:>8.2f} | {m_link:>9.2f} | {pct:>7.1f}%")
        except Exception as e:
            print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    run_experiment()
