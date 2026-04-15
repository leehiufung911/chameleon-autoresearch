
import sys
import os
import math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.metrics import roc_auc_score

# Add chameleon_local to path
sys.path.append("chameleon_local")
import chameleon

def run_experiment():
    print("Loading data...", flush=True)
    df_labels = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    df_user["label"] = "user"
    
    df_full = pd.concat([df_labels, df_user], ignore_index=True).sample(20, random_state=42)
    
    n_conf = 50
    results = []

    print(f"Processing {len(df_full)} molecules with N_CONF={n_conf}...", flush=True)

    for i, row in df_full.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        print(f"  [{i+1}/{len(df_full)}] {name}...", end="", flush=True)
        
        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None:
            print(" (Failed)")
            continue
            
        try:
            mol, cids = chameleon.embed_conformers(mol0, n_conf=n_conf, n_threads=0)
            energies = chameleon.minimize(mol, n_threads=0)
            keep, k_energies = chameleon.butina_prune(mol, energies, rms_cut=0.75)
            
            psa_list = []
            rg_list = []
            imhb_list = []
            radii = chameleon.rdFreeSASA.classifyAtoms(mol)
            topo = chameleon.bond_path_distance_matrix(mol)
            
            for cid in keep:
                psa_list.append(chameleon.compute_3d_psa(mol, cid, radii))
                rg_list.append(chameleon.compute_rg(mol, cid))
                imhb_list.append(chameleon.compute_imhb(mol, cid, topo))
                
            psa = np.array(psa_list)
            rg = np.array(rg_list)
            imhb = np.array(imhb_list)
            
            mw = Descriptors.MolWt(mol0)
            psa_range = psa.max() - psa.min()
            psa_ratio = psa.max() / max(psa.min(), 1.0)
            rg_ratio = rg.max() / max(rg.min(), 0.1)
            imhb_mean = imhb.mean()
            
            # 1. Original CI
            term_psa = psa_range / math.sqrt(mw)
            term_rg = math.log(rg_ratio)
            ci_orig = 2.0 * term_psa + 1.0 * imhb_mean + 3.0 * term_rg
            
            # 2. Ratio CI (Scaled to be comparable)
            # Typical psa_ratio is 1.1 - 1.5. (Ratio - 1.0) * 10 is ~1.0 - 5.0.
            term_psa_ratio = (psa_ratio - 1.0) * 10.0
            ci_ratio = 1.0 * term_psa_ratio + 1.0 * imhb_mean + 3.0 * term_rg

            results.append({
                "name": name,
                "label": label,
                "is_pos": (label in ("chameleon", "protac")) or (name in ("protac_1", "protac_2")),
                "ci_orig": ci_orig,
                "ci_ratio": ci_ratio,
                "psa_ratio": psa_ratio,
                "imhb_mean": imhb_mean,
                "term_psa": term_psa,
                "term_rg": term_rg
            })
            print(" Done.")
            
        except Exception as e:
            print(f" (Error: {e})")
            continue

    df_res = pd.DataFrame(results)
    
    # Evaluation on benchmark
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    print("\n--- Benchmark AUC (N=49) ---")
    y = df_eval["is_pos"].astype(int)
    for col in ["ci_orig", "ci_ratio", "psa_ratio", "imhb_mean", "term_psa"]:
        auc = roc_auc_score(y, df_eval[col])
        print(f"{col:15}: AUC = {auc:.3f}")

    print("\n--- User PROTACs (Verdict at threshold 4.0) ---")
    df_u = df_res[df_res["label"] == "user"].copy()
    print(df_u[["name", "ci_orig", "ci_ratio", "imhb_mean"]])

if __name__ == "__main__":
    run_experiment()
