
import sys
import os
import math
import numpy as np
from rdkit import Chem
sys.path.append(os.getcwd())
from chameleon_local.chameleon import summarize, Summary

def get_ring_systems(mol):
    ring_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()]
    if not ring_atoms: return []
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
    queue = []
    visited = {}
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
    if found_target is None: return []
    path = []
    curr = found_target
    while curr is not None:
        path.append(curr)
        curr = visited[curr]
    return path[::-1]

def get_linker_stats(mol):
    systems = get_ring_systems(mol)
    if len(systems) < 2: return 0, 0, 0.0
    max_non_ring = -1
    best_path_stats = (0, 0, 0.0)
    for i in range(len(systems)):
        for j in range(i + 1, len(systems)):
            path = get_shortest_path_between_systems(mol, systems[i], systems[j])
            if path:
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
    test_mols = [
        ("protac_1", "O=C(C(N1C(C2=CC=CC(NCCOCCOCCOCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=C2C1=O)=O)CC7)NC7=O", "chameleon"),
        ("protac_2", "O=C(C(N1C(C2=CC=CC(NC(COCCOCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O", "chameleon"),
        ("protac_3", "O=C(C(N1C(C2=CC=CC(OCC(NCCCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O", "nonchameleon"),
        ("dBET6", "CC1=C(SC2=C1C(=NC(C3=NN=C(N32)C)CC(=O)NCCCCCCCCNC(=O)COC4=CC=CC5=C4C(=O)N(C5=O)C6CCC(=O)NC6=O)C7=CC=C(C=C7)Cl)C", "chameleon"),
        ("MZ1", "Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC(=O)COCCOCCOCCNC(=O)C[C@@H]2N=C(c3ccc(Cl)cc3)c3c(sc(C)c3C)-n3c(C)nnc32)C(C)(C)C)cc1", "chameleon"),
        ("cyclosporinA", "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C", "chameleon"),
        ("atorvastatin", "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4", "nonchameleon")
    ]
    
    print(f"{'Name':<15} | {'GT':<10} | {'Density':>7} | {'IMHB':>6} | {'CI_orig':>7} | {'CI_bal':>7} | {'Correct?'}")
    print("-" * 85)
    
    for name, smi, label in test_mols:
        mol = Chem.MolFromSmiles(smi)
        n_total, n_hetero, density = get_linker_stats(mol)
        s = summarize(name, smi, n_conf=150, verbose=False)
        
        ci_orig = s.chameleonic_index
        imhb = s.imhb_mean
        
        # Balance score logic
        ci_bal = ci_orig
        penalty_reason = "None"
        if density > 0 and density < 0.25 and imhb < 0.9:
            ci_bal = ci_orig * 0.5
            penalty_reason = "AlkylLinker"
        
        # Verdict
        verdict = "CHAM" if ci_bal >= 2.0 else "NON"
        gt_short = "CHAM" if label == "chameleon" else "NON"
        correct = "YES" if verdict == gt_short else "NO"
        
        print(f"{name:<15} | {gt_short:<10} | {density:7.3f} | {imhb:6.2f} | {ci_orig:7.2f} | {ci_bal:7.2f} | {correct} ({penalty_reason})")

if __name__ == "__main__":
    run_experiment()
