import sys
import os
import math
import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

sys.path.append("chameleon_local")
import chameleon

def run_experiment():
    print("Loading data...", flush=True)
    df_labels = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    df_user["label"] = "user"
    
    chams = df_labels[df_labels["label"] == "chameleon"].sample(10, random_state=42)
    nons = df_labels[df_labels["label"] == "nonchameleon"].sample(10, random_state=42)
    
    df_test = pd.concat([chams, nons, df_user], ignore_index=True)
    n_conf = 30
    results = []

    print(f"Processing {len(df_test)} molecules with N_CONF={n_conf}...", flush=True)

    for i, row in df_test.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        print(f"[{i+1}/{len(df_test)}] {name}...", end=" ", flush=True)
        
        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None:
            print("(Failed)")
            continue
            
        try:
            mol, cids = chameleon.embed_conformers(mol0, n_conf=n_conf, n_threads=0)
            energies = chameleon.minimize(mol, n_threads=0)
            keep, _ = chameleon.butina_prune(mol, energies, rms_cut=0.5)
            
            if len(keep) < 3:
                print(f"(Too few confs: {len(keep)})")
                continue

            psa_list = []
            rg_list = []
            radii = chameleon.rdFreeSASA.classifyAtoms(mol)
            
            for cid in keep:
                psa_list.append(chameleon.compute_3d_psa(mol, cid, radii))
                rg_list.append(chameleon.compute_rg(mol, cid))
                
            psa = np.array(psa_list)
            rg = np.array(rg_list)
            
            corr_rg_psa, _ = spearmanr(rg, psa)
            if np.isnan(corr_rg_psa):
                corr_rg_psa = 0.0

            is_pos = (label in ("chameleon", "protac")) or (name in ("protac_1", "protac_2"))
            
            results.append({
                "name": name,
                "label": label,
                "is_pos": is_pos,
                "corr_rg_psa": corr_rg_psa,
            })
            print(f"corr={corr_rg_psa:.3f}")
            
        except Exception as e:
            print(f"(Error: {e})")
            continue

    df_res = pd.DataFrame(results)
    
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    print("\n--- Benchmark AUC ---")
    y = df_eval["is_pos"].astype(int)
    auc = roc_auc_score(y, df_eval["corr_rg_psa"])
    print(f"corr_rg_psa AUC: {auc:.3f}")

    print("\n--- User PROTACs ---")
    for _, r in df_res[df_res["label"] == "user"].iterrows():
        print(f"{r['name']:15}: corr={r['corr_rg_psa']:.3f} (True: {r['is_pos']})")

if __name__ == "__main__":
    run_experiment()
