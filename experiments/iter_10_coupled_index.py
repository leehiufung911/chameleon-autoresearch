
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import csv

# Add chameleon_local to path
sys.path.append("chameleon_local")
import chameleon

def roc_auc(scores, positives):
    """Simple AUC implementation."""
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
    if n_pos == 0 or n_neg == 0:
        return 0.5
    sum_ranks_pos = sum(r for r, p in zip(ranks, positives) if p)
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)

def run_experiment():
    # Load labels
    labels = {}
    with open("chameleon_local/labelled_set.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            labels[row["name"]] = row["label"]

    # Fixed small balanced subset
    all_label_data = {}
    with open("chameleon_local/labelled_set.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            all_label_data[row["name"]] = row["smiles"]

    pos_names = ["tacrolimus", "sirolimus", "everolimus", "cyclosporinA", "rifampicin"]
    neg_names = ["aspirin", "ibuprofen", "acetaminophen", "caffeine", "diazepam", "warfarin", "atenolol"]
    
    all_items = []
    for n in pos_names + neg_names:
        if n in all_label_data:
            all_items.append((n, all_label_data[n]))

    with open("chameleon_local/user_protacs.tsv", "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) == 2: all_items.append((row[0], row[1]))

    N_CONF = 15
    results = []
    
    print(f"{'Name':<20} | {'IMHB':>5} | {'Nrot':>4} | {'HBD':>3} | {'Index':>8} | {'Status'}")
    print("-" * 65)

    for name, smi in all_items:
        try:
            summary = chameleon.summarize(name, smi, n_conf=N_CONF, n_threads=1, verbose=False)
            mol = Chem.MolFromSmiles(smi)
            n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            if hbd == 0: hbd = 0.5
            
            # IMHB-Entropy-Saturation Coupled Index
            index = summary.imhb_mean * math.log(2 + n_rot) / hbd
            
            label = labels.get(name)
            is_pos = (label in ("chameleon", "protac")) or (name in ("protac_1", "protac_2"))
            
            results.append({
                "name": name,
                "is_pos": is_pos,
                "index": index,
                "ci": summary.chameleonic_index
            })
            
            status = "POS" if is_pos else "NEG"
            if name.startswith("protac"):
                print(f"**{name[:18]:<18} | {summary.imhb_mean:>5.2f} | {int(n_rot):>4d} | {int(hbd):>3d} | {index:>8.4f} | {status}**")
            else:
                print(f"{name[:20]:<20} | {summary.imhb_mean:>5.2f} | {int(n_rot):>4d} | {int(hbd):>3d} | {index:>8.4f} | {status}")
                
        except Exception as e:
            pass

    # Calculate AUCs
    ys = [r["is_pos"] for r in results]
    scores_ci = [r["ci"] for r in results]
    scores_idx = [r["index"] for r in results]

    auc_ci = roc_auc(scores_ci, ys)
    auc_idx = roc_auc(scores_idx, ys)

    print("\n" + "="*50)
    print(f"Original CI AUC:       {auc_ci:.3f}")
    print(f"Coupled Index AUC:     {auc_idx:.3f}")
    print("="*50)

if __name__ == "__main__":
    run_experiment()
