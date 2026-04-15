
import sys
import os
import math
from rdkit import Chem

def get_ring_systems(mol):
    """Find connected components of ring atoms."""
    ring_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()]
    if not ring_atoms:
        return []
    
    systems = []
    visited = set()
    for idx in ring_atoms:
        if idx not in visited:
            comp = []
            stack = [idx]
            visited.add(idx)
            while stack:
                curr = stack.pop()
                comp.append(curr)
                for neighbor in mol.GetAtomWithIdx(curr).GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if neighbor.IsInRing() and n_idx not in visited:
                        visited.add(n_idx)
                        stack.append(n_idx)
            systems.append(set(comp))
    return systems

def get_shortest_path_between_systems(mol, set1, set2):
    """Find shortest path using BFS to avoid GetDistanceMatrix/numpy."""
    # BFS starting from all atoms in set1
    queue = []
    visited = {} # atom_idx -> parent_idx
    for idx in set1:
        queue.append(idx)
        visited[idx] = None
    
    found_target = None
    curr_idx_in_queue = 0
    while curr_idx_in_queue < len(queue):
        curr = queue[curr_idx_in_queue]
        curr_idx_in_queue += 1
        
        if curr in set2:
            found_target = curr
            break
            
        for neighbor in mol.GetAtomWithIdx(curr).GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited:
                visited[n_idx] = curr
                queue.append(n_idx)
    
    if found_target is None:
        return []
    
    # Reconstruct path
    path = []
    curr = found_target
    while curr is not None:
        path.append(curr)
        curr = visited[curr]
    return path[::-1]

def get_linker_stats(mol):
    """Find the path between ring systems that has the most non-ring atoms."""
    systems = get_ring_systems(mol)
    if len(systems) < 2:
        return 0, 0, 0.0 # No linker between systems
    
    max_non_ring = -1
    best_path_stats = (0, 0, 0.0) # n_total, n_hetero, density
    
    for i in range(len(systems)):
        for j in range(i + 1, len(systems)):
            path = get_shortest_path_between_systems(mol, systems[i], systems[j])
            
            if path:
                # Linker atoms are those in the path excluding the endpoints (which are in rings)
                # Actually, the endpoints of the shortest path might be in rings.
                # We want the non-ring segment.
                non_ring_segment = [idx for idx in path if not mol.GetAtomWithIdx(idx).IsInRing()]
                n_total = len(non_ring_segment)
                if n_total > 0:
                    n_hetero = 0
                    for idx in non_ring_segment:
                        atom = mol.GetAtomWithIdx(idx)
                        if atom.GetSymbol() in ("N", "O", "S", "P", "F", "Cl"):
                            n_hetero += 1
                    density = n_hetero / n_total
                    
                    if n_total > max_non_ring:
                        max_non_ring = n_total
                        best_path_stats = (n_total, n_hetero, density)
                        
    return best_path_stats

def run_experiment():
    # Load user PROTACs
    user_protacs = []
    try:
        with open("chameleon_local/user_protacs.tsv") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split("\t")
                if len(parts) == 2:
                    user_protacs.append((parts[0], parts[1]))
    except Exception as e:
        print(f"Error loading user PROTACs: {e}")
            
    # Load benchmark
    benchmark = []
    try:
        with open("chameleon_local/labelled_set.tsv") as f:
            header = f.readline()
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split("\t")
                if len(parts) == 3:
                    benchmark.append((parts[0], parts[1], parts[2]))
    except Exception as e:
        print(f"Error loading benchmark: {e}")
            
    print(f"{'Name':<20} | {'Label':<12} | {'Len':<4} | {'Het':<4} | {'Density':<7}")
    print("-" * 60)
    
    results = []
    
    for name, smi in user_protacs:
        mol = Chem.MolFromSmiles(smi)
        if mol is None: continue
        n_total, n_hetero, density = get_linker_stats(mol)
        label = "user_cham" if "3" not in name else "user_non"
        print(f"{name:<20} | {label:<12} | {n_total:<4} | {n_hetero:<4} | {density:7.3f}")
        results.append((name, label, density))

    print("-" * 60)
    
    for label, name, smi in benchmark:
        mol = Chem.MolFromSmiles(smi)
        if mol is None: continue
        n_total, n_hetero, density = get_linker_stats(mol)
        # Only show PROTACs and a few interesting ones
        if label == "protac" or density > 0.3: # Filter for interesting ones
            print(f"{name:<20} | {label:<12} | {n_total:<4} | {n_hetero:<4} | {density:7.3f}")
            results.append((name, label, density, label))

    # Analysis
    print("\nSummary:")
    p1 = [r for r in results if r[0] == "protac_1"]
    p3 = [r for r in results if r[0] == "protac_3"]
    if p1 and p3:
        print(f"protac_1 (Cham) Density: {p1[0][2]:.3f}")
        print(f"protac_3 (Non)  Density: {p3[0][2]:.3f}")
    
    dbet6 = [r for r in results if r[0] == "dBET6"]
    if dbet6:
        print(f"dBET6 (Cham) Density: {dbet6[0][2]:.3f}")

if __name__ == "__main__":
    run_experiment()
