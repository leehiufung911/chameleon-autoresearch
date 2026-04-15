
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from sklearn.metrics import roc_auc_score

def count_ethers(mol):
    # SMARTS for aliphatic ether oxygen
    # O, 0 H, valence 2, not double bonded to C, N, S
    ether_smarts = "[O;H0;v2;!$(O=C);!$(O=P);!$(O=N);!$(O=S)]"
    patt = Chem.MolFromSmarts(ether_smarts)
    return len(mol.GetSubstructMatches(patt))

def run():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    results = []
    print("Processing User PROTACs...", flush=True)
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        n_ethers = count_ethers(mol)
        results.append({
            "name": name, "label": "protac", "n_ethers": n_ethers
        })
        print(f"  {name}: n_ethers={n_ethers}", flush=True)

    print("Processing Benchmark...", flush=True)
    for i, row in df_bench.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        mol = Chem.MolFromSmiles(smiles)
        n_ethers = count_ethers(mol)
        results.append({
            "name": name, "label": label, "n_ethers": n_ethers
        })

    df_res = pd.DataFrame(results)
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
    
    # AUC for ethers (expect positive correlation with chameleonicity in PROTACs)
    # But wait, cyclosporin has no ethers!
    auc = roc_auc_score(df_eval["target"], df_eval["n_ethers"])
    print("\n--- AUC RESULTS (Ether Count) ---", flush=True)
    print("AUC: {:.3f}".format(auc))
    
    print("\n--- Summary Stats ---", flush=True)
    print(df_eval.groupby("label")["n_ethers"].mean())
    
    print("\n--- User PROTACs Detail ---")
    print(df_res[df_res["name"].str.startswith("protac")][["name", "n_ethers"]])

if __name__ == "__main__":
    run()
