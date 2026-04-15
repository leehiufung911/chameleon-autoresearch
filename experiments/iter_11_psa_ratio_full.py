import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from sklearn.metrics import roc_auc_score

# Add chameleon_local to path
sys.path.append("chameleon_local")
import chameleon

def run_experiment():
    print("Loading data...", flush=True)
    # Benchmark
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    # User PROTACs
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    df_user["label"] = "user"

    df_full = pd.concat([df_bench, df_user], ignore_index=True)

    n_conf = 10
    results = []

    print(f"Starting processing of {len(df_full)} molecules with n_conf={n_conf}...", flush=True)

    for i, row in df_full.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        print(f"  {i+1}/{len(df_full)}: {name}", flush=True)
        
        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None:
            continue
            
        # Conformer generation & MMFF minimization
        try:
            mol, cids = chameleon.embed_conformers(mol0, n_conf=n_conf, n_threads=0)
            energies = chameleon.minimize(mol, n_threads=0)
            keep, _ = chameleon.butina_prune(mol, energies)
        except Exception as e:
            print(f"Error processing {name}: {e}")
            continue

        psa_list = []
        radii = chameleon.rdFreeSASA.classifyAtoms(mol)
        
        for cid in keep:
            psa_list.append(chameleon.compute_3d_psa(mol, cid, radii))
            
        psa = np.array(psa_list)
        
        # Raw descriptors
        psa_max = psa.max()
        psa_min = psa.min()
        psa_ratio = psa_max / max(psa_min, 1.0)
        
        results.append({
            "name": name,
            "label": label,
            "psa_max": psa_max,
            "psa_min": psa_min,
            "psa_ratio": psa_ratio
        })

    df_res = pd.DataFrame(results)
    
    print("\n--- Benchmark Set Evaluation (N=49) ---")
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    if len(df_eval) > 0:
        df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
        auc = roc_auc_score(df_eval["target"], df_eval["psa_ratio"])
        print(f"PSA_ratio AUC: {auc:.3f}")
    else:
        print("Not enough data for AUC")

    print("\n--- User PROTACs ---")
    df_user_res = df_res[df_res["label"] == "user"].copy()
    print(df_user_res[["name", "psa_max", "psa_min", "psa_ratio"]])

if __name__ == "__main__":
    try:
        run_experiment()
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)
