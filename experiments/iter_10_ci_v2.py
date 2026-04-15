
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))
import chameleon

def count_potential_imhb_pairs(mol):
    topo = Chem.GetDistanceMatrix(mol)
    donors = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in ("N", "O"):
            for n in a.GetNeighbors():
                if n.GetSymbol() == "H":
                    donors.append(a.GetIdx())
                    break
    acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ("N", "O")]
    
    pairs = 0
    for d in donors:
        for a in acceptors:
            if d != a and topo[d, a] >= 4:
                pairs += 1
    return max(pairs, 1)

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

def roc_auc(scores, positives):
    n = len(scores)
    order = sorted(range(n), key=lambda i: scores[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    n_pos = sum(1 for p in positives if p)
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0: return 0.5
    sum_ranks_pos = sum(r for r, p in zip(ranks, positives) if p)
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)

def run_experiment():
    labels = {}
    with open("chameleon_local/labelled_set.tsv") as f:
        header = f.readline().strip().split("\t")
        iL = header.index("label")
        iN = header.index("name")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                labels[parts[iN]] = parts[iL]

    all_items = []
    with open("chameleon_local/benchmark.tsv") as f:
        for i, line in enumerate(f):
            if i > 30: break
            parts = line.strip().split("\t")
            if len(parts) == 2:
                all_items.append((parts[0], parts[1]))
    
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                all_items.append((parts[0], parts[1]))

    N_CONF = 10
    print(f"{'Name':<20} | {'CI_orig':>7} | {'Sat':>7} | {'MaxC':>4} | {'CI_v2':>7} | {'Label'}")
    print("-" * 70)

    scores_orig = []
    scores_v2 = []
    ys = []

    for name, smi in all_items:
        try:
            mol = Chem.MolFromSmiles(smi)
            mol_hs = Chem.AddHs(mol)
            
            # 1. Run standard pipeline
            summary = chameleon.summarize(name, smi, n_conf=N_CONF, verbose=False)
            ci_orig = summary.chameleonic_index
            
            # 2. Compute Saturation
            pot = count_potential_imhb_pairs(mol_hs)
            sat = summary.imhb_mean / pot
            
            # 3. Compute Max Carbon Path
            max_c = get_max_nonring_carbon_path(mol)
            
            # CI v3: Additive CI + Gated MaxC Penalty
            ci_v3 = ci_orig
            
            # Gated Aliphatic Penalty
            if max_c >= 5 and summary.imhb_mean < 1.0:
                ci_v3 *= 0.5
            
            label = labels.get(name)
            is_pos = (label in ("chameleon", "protac")) or (name in ("protac_1", "protac_2"))
            
            if label or name.startswith("protac"):
                ys.append(is_pos)
                scores_orig.append(ci_orig)
                scores_v2.append(ci_v3) # Rename v3 to v2 for consistency with plot
            
            l_str = "CHAM" if is_pos else "NON"
            print(f"{name[:20]:<20} | {ci_orig:>7.2f} | {max_c:>4} | {summary.imhb_mean:>5.2f} | {ci_v3:>7.2f} | {l_str}")
            
        except Exception as e:
            # print(f"Error processing {name}: {e}")
            pass

    auc_orig = roc_auc(scores_orig, ys)
    auc_v2 = roc_auc(scores_v2, ys)

    print("\n" + "="*50)
    print(f"Original AUC (N={N_CONF}): {auc_orig:.3f}")
    print(f"CI v2 AUC (Sat+Veto): {auc_v2:.3f}")
    print("="*50)

if __name__ == "__main__":
    run_experiment()
