
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
from sklearn.metrics import roc_auc_score

sys.path.append("chameleon_local")
import chameleon

def count_hbd(mol):
    return rdMolDescriptors.CalcNumHBD(mol)

def run():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    N_CONF = 30
    results = []
    
    print("Processing User PROTACs...", flush=True)
    for name, smiles in zip(df_user["name"], df_user["smiles"]):
        mol0 = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol0)
        AllChem.EmbedMultipleConfs(mol, numConfs=N_CONF, randomSeed=42)
        topo = chameleon.bond_path_distance_matrix(mol)
        imhb_pre_mean = np.mean([chameleon.compute_imhb(mol, cid, topo) for cid in range(mol.GetNumConformers())])
        s = chameleon.summarize(name, smiles, n_conf=N_CONF, verbose=False)
        hbd = count_hbd(mol0)
        results.append({
            "name": name, "label": "protac", "hbd": hbd, 
            "imhb_pre": imhb_pre_mean,
            "sat_ratio_pre": imhb_pre_mean / hbd if hbd > 0 else 1.0,
            "ci": s.chameleonic_index
        })
        print(f"  {name}: SatRatio_pre={results[-1]['sat_ratio_pre']:.3f}, CI={s.chameleonic_index:.2f}", flush=True)

    print("\nProcessing Benchmark (subset)...", flush=True)
    df_sub = pd.concat([df_bench[df_bench["label"] == "chameleon"].head(5), 
                        df_bench[df_bench["label"] == "nonchameleon"].head(5)])
    
    for name, smiles, label in zip(df_sub["name"], df_sub["smiles"], df_sub["label"]):
        mol0 = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol0)
        AllChem.EmbedMultipleConfs(mol, numConfs=N_CONF, randomSeed=42)
        topo = chameleon.bond_path_distance_matrix(mol)
        imhb_pre_mean = np.mean([chameleon.compute_imhb(mol, cid, topo) for cid in range(mol.GetNumConformers())])
        s = chameleon.summarize(name, smiles, n_conf=N_CONF, verbose=False)
        hbd = count_hbd(mol0)
        results.append({
            "name": name, "label": label, "hbd": hbd, 
            "imhb_pre": imhb_pre_mean,
            "sat_ratio_pre": imhb_pre_mean / hbd if hbd > 0 else 1.0,
            "ci": s.chameleonic_index
        })
        print(f"  {name}: SatRatio_pre={results[-1]['sat_ratio_pre']:.3f}", flush=True)

    df_res = pd.DataFrame(results)
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
    
    print("\n--- AUC RESULTS (Partial Benchmark, N_CONF=30) ---", flush=True)
    print("AUC CI:            {:.3f}".format(roc_auc_score(df_eval["target"], df_eval["ci"])))
    print("AUC Sat_Ratio_pre: {:.3f}".format(roc_auc_score(df_eval["target"], df_eval["sat_ratio_pre"])))
    
    print("\n--- User PROTACs Detail ---")
    print(df_res[df_res["label"] == "protac"][["name", "hbd", "imhb_pre", "sat_ratio_pre", "ci"]])

if __name__ == "__main__":
    run()
