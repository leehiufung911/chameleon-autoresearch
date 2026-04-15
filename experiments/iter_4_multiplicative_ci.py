
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from sklearn.metrics import roc_auc_score

sys.path.append("chameleon_local")
import chameleon

def run():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    results = []
    
    # Process User PROTACs
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        s = chameleon.summarize(name, smiles, n_conf=5, verbose=False)
        psa_range = s.psa3d_max - s.psa3d_min
        mw = s.mw
        rg_ratio = s.rg_max / max(s.rg_min, 1e-6)
        term_psa = psa_range / np.sqrt(mw)
        term_imhb = s.imhb_mean
        term_rg = np.log(max(rg_ratio, 1.0))
        ci_orig = 2.0 * term_psa + 1.0 * term_imhb + 3.0 * term_rg
        ci_mult = (term_psa + 3.0 * term_rg) * term_imhb
        results.append({
            "name": name, "label": "protac", "ci_orig": ci_orig, "ci_mult": ci_mult, "imhb": term_imhb
        })

    # Process 10 Benchmark
    df_eval_set = pd.concat([
        df_bench[df_bench["label"] == "chameleon"].sample(5, random_state=42),
        df_bench[df_bench["label"] == "nonchameleon"].sample(5, random_state=42)
    ])
    for i, row in df_eval_set.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        try:
            s = chameleon.summarize(name, smiles, n_conf=5, verbose=False)
            psa_range = s.psa3d_max - s.psa3d_min
            mw = s.mw
            rg_ratio = s.rg_max / max(s.rg_min, 1e-6)
            term_psa = psa_range / np.sqrt(mw)
            term_imhb = s.imhb_mean
            term_rg = np.log(max(rg_ratio, 1.0))
            ci_orig = 2.0 * term_psa + 1.0 * term_imhb + 3.0 * term_rg
            ci_mult = (term_psa + 3.0 * term_rg) * term_imhb
            results.append({
                "name": name, "label": label, "ci_orig": ci_orig, "ci_mult": ci_mult, "imhb": term_imhb
            })
            print(f"  {name} ok", flush=True)
        except: pass

    df_res = pd.DataFrame(results)
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
    
    auc_orig = roc_auc_score(df_eval["target"], df_eval["ci_orig"])
    auc_mult = roc_auc_score(df_eval["target"], df_eval["ci_mult"])
    
    print("\n--- AUC RESULTS (N_CONF=5, sample=10) ---")
    print(f"AUC CI_orig: {auc_orig:.3f}")
    print(f"AUC CI_mult: {auc_mult:.3f}")
    
    print("\n--- User PROTACs Detail ---")
    print(df_res[df_res["label"] == "protac"][["name", "ci_orig", "ci_mult", "imhb"]])

if __name__ == "__main__":
    run()
