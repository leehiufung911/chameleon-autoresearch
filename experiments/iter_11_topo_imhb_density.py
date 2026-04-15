import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from sklearn.metrics import roc_auc_score

sys.path.append("chameleon_local")
import chameleon

def run_experiment():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    df_user["label"] = "user"

    df_full = pd.concat([df_bench, df_user], ignore_index=True)
    
    n_conf = 10
    results = []

    print(f"Processing {len(df_full)} molecules (n_conf={n_conf})...", flush=True)

    for i, row in df_full.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        print(f"  {i+1}/{len(df_full)}: {name}", flush=True)
        
        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None: continue
            
        try:
            mol, cids = chameleon.embed_conformers(mol0, n_conf=n_conf, n_threads=0)
            energies = chameleon.minimize(mol, n_threads=0)
            keep, _ = chameleon.butina_prune(mol, energies)
        except Exception as e:
            print(f"Error on {name}: {e}")
            continue

        topo = chameleon.bond_path_distance_matrix(mol)
        
        imhb_all = []
        imhb_long = []
        
        for cid in keep:
            imhb_all.append(chameleon.compute_imhb(mol, cid, topo, min_topo=4))
            imhb_long.append(chameleon.compute_imhb(mol, cid, topo, min_topo=11))
            
        imhb_all = np.array(imhb_all)
        imhb_long = np.array(imhb_long)
        
        if len(imhb_all) == 0:
            continue
            
        imhb_mean_all = imhb_all.mean()
        imhb_mean_long = imhb_long.mean()
        
        results.append({
            "name": name,
            "label": label,
            "imhb_mean_all": imhb_mean_all,
            "imhb_mean_long": imhb_mean_long
        })

    df_res = pd.DataFrame(results)
    
    print("\n--- Benchmark Set Evaluation (N=49) ---")
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    if len(df_eval) > 0:
        df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
        auc_all = roc_auc_score(df_eval["target"], df_eval["imhb_mean_all"])
        auc_long = roc_auc_score(df_eval["target"], df_eval["imhb_mean_long"])
        print(f"IMHB (min_topo=4)  AUC: {auc_all:.3f}")
        print(f"IMHB (min_topo=11) AUC: {auc_long:.3f}")
    
    print("\n--- User PROTACs ---")
    df_user_res = df_res[df_res["label"] == "user"].copy()
    print(df_user_res[["name", "imhb_mean_all", "imhb_mean_long"]])

if __name__ == "__main__":
    run_experiment()
