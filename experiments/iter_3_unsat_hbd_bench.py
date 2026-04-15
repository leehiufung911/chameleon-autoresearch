
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sklearn.metrics import roc_auc_score

sys.path.append("chameleon_local")
import chameleon

def run():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    # Sample 5 of each label for speed
    df_cham = df_bench[df_bench["label"] == "chameleon"].sample(5, random_state=42)
    df_non = df_bench[df_bench["label"] == "nonchameleon"].sample(5, random_state=42)
    df_eval_set = pd.concat([df_cham, df_non])
    
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    N_CONF = 5
    results = []
    
    print(f"Processing {len(df_user)} User PROTACs...", flush=True)
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        try:
            s = chameleon.summarize(name, smiles, n_conf=N_CONF, verbose=False)
            hbd = rdMolDescriptors.CalcNumHBD(Chem.MolFromSmiles(smiles))
            results.append({
                "name": name, "label": "protac", "hbd": hbd, "imhb": s.imhb_mean, 
                "unsat": hbd - s.imhb_mean, "ci": s.chameleonic_index
            })
            print(f"  {name}: Unsat={results[-1]['unsat']:.2f}, CI={results[-1]['ci']:.2f}", flush=True)
        except: pass

    print(f"Processing {len(df_eval_set)} Benchmark molecules...", flush=True)
    for i, row in df_eval_set.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        try:
            s = chameleon.summarize(name, smiles, n_conf=N_CONF, verbose=False)
            hbd = rdMolDescriptors.CalcNumHBD(Chem.MolFromSmiles(smiles))
            results.append({
                "name": name, "label": label, "hbd": hbd, "imhb": s.imhb_mean, 
                "unsat": hbd - s.imhb_mean, "ci": s.chameleonic_index
            })
            if (i+1) % 5 == 0:
                print(f"  Progress: {i+1}/{len(df_eval_set)}", flush=True)
        except: pass

    df_res = pd.DataFrame(results)
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
    
    print("\n--- AUC (N_CONF=10, sample=30) ---", flush=True)
    print("AUC CI (baseline): {:.3f}".format(roc_auc_score(df_eval["target"], df_eval["ci"])))
    print("AUC Unsat_HBD:    {:.3f}".format(roc_auc_score(df_eval["target"], -df_eval["unsat"])))
    
    print("\n--- Summary Stats ---", flush=True)
    print(df_eval.groupby("label")[["unsat", "ci"]].mean())

if __name__ == "__main__":
    run()
