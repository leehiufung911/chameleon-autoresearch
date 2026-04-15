
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sklearn.metrics import roc_auc_score

def get_topological_imhb_potential(mol, min_dist=7, max_dist=15):
    # Donors: H-bond donors (N-H, O-H)
    # Acceptors: H-bond acceptors (N, O)
    
    # Get atom indices for donors and acceptors
    donors = []
    acceptors = []
    
    # We'll use simple definitions for this experiment
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ["N", "O"]:
            # Acceptor
            acceptors.append(atom.GetIdx())
            # Donor if it has at least one H
            if atom.GetTotalNumHs() > 0:
                donors.append(atom.GetIdx())
    
    if not donors or not acceptors:
        return 0
    
    # Get topological distance matrix
    dist_mat = Chem.GetDistanceMatrix(mol)
    
    count = 0
    for d in donors:
        for a in acceptors:
            if d == a: continue
            dist = dist_mat[d, a]
            if min_dist <= dist <= max_dist:
                count += 1
    return count

def run():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    # Load original CI scores for benchmark if available
    # Actually I'll just use my 2D score
    
    results = []
    print("Processing User PROTACs...", flush=True)
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        topo_score = get_topological_imhb_potential(mol)
        results.append({
            "name": name, "label": "protac", "topo_score": topo_score
        })
        print(f"  {name}: topo_score={topo_score}", flush=True)

    print("Processing Benchmark...", flush=True)
    for i, row in df_bench.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        mol = Chem.MolFromSmiles(smiles)
        topo_score = get_topological_imhb_potential(mol)
        results.append({
            "name": name, "label": label, "topo_score": topo_score
        })

    df_res = pd.DataFrame(results)
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
    
    auc = roc_auc_score(df_eval["target"], df_eval["topo_score"])
    print("\n--- AUC RESULTS (Topological IMHB Potential) ---", flush=True)
    print("AUC: {:.3f}".format(auc))
    
    print("\n--- Summary Stats ---", flush=True)
    print(df_eval.groupby("label")["topo_score"].mean())
    
    print("\n--- User PROTACs Detail ---")
    print(df_res[df_res["label"] == "protac"][["name", "topo_score"]])

if __name__ == "__main__":
    run()
