
import os
import sys

# Manually add conda env paths for Windows
env_path = r"C:\Users\mic23\miniconda3\envs\chameleon"
os.environ["PATH"] = (
    os.path.join(env_path, "Library", "bin") + ";" +
    os.path.join(env_path, "Scripts") + ";" +
    env_path + ";" +
    os.environ.get("PATH", "")
)

import pandas as pd
import numpy as np
import math
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors

# Add current dir to path to import chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))
import chameleon as cham

def get_sat_ratio_pre(mol, n_conf=30):
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    params.numThreads = 1
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    if not cids: return 0.0
    
    total_hbd = rdMolDescriptors.CalcNumHBD(mol)
    if total_hbd <= 0: return 1.0 # Satisfied or N/A
    
    topo = cham.bond_path_distance_matrix(mol)
    imhb_counts = []
    for cid in cids:
        imhb_counts.append(cham.compute_imhb(mol, cid, topo))
    
    mean_imhb = np.mean(imhb_counts)
    return mean_imhb / total_hbd

def run_experiment():
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    df_user["label"] = "protac"
    user_truth = {"protac_1": "chameleon", "protac_2": "chameleon", "protac_3": "nonchameleon"}

    results = []
    bench_to_run = [
        "cyclosporinA", "tacrolimus", "sirolimus", "erythromycin", "clarithromycin", 
        "aspirin", "ibuprofen", "atenolol", "metoprolol", "propranolol", "nifedipine"
    ]
    
    print("Processing molecules...", flush=True)
    for df in [df_bench, df_user]:
        for i, row in df.iterrows():
            name = row["name"]
            if df is df_bench and name not in bench_to_run: continue
            
            smiles = row["smiles"]
            true_label = row["label"]
            if true_label == "protac":
                true_label = user_truth.get(name, "unknown")
            
            print(f"  {name}...", end=" ", flush=True)
            try:
                s = cham.summarize(name, smiles, n_conf=30, verbose=False)
                mol = Chem.MolFromSmiles(smiles)
                sat_ratio = get_sat_ratio_pre(mol, n_conf=30)
                
                results.append({
                    "name": name,
                    "label": true_label,
                    "CI_orig": s.chameleonic_index,
                    "SatRatio_pre": sat_ratio
                })
                print(f"CI={s.chameleonic_index:.2f}, SatRatio={sat_ratio:.3f}")
            except Exception as e:
                print(f"FAILED: {e}")
                
    res_df = pd.DataFrame(results)
    
    # Final proposed CI: Veto if SatRatio_pre is extremely low
    res_df["CI_final"] = res_df["CI_orig"] * np.where(res_df["SatRatio_pre"] > 0.005, 1.0, 0.4)
    
    from sklearn.metrics import roc_auc_score
    bench_only = res_df[res_df["label"] != "unknown"].copy()
    bench_only["y_true"] = bench_only["label"].map({"chameleon": 1, "nonchameleon": 0})
    bench_only = bench_only.dropna(subset=["y_true"])
    
    if len(bench_only["y_true"].unique()) > 1:
        auc_orig = roc_auc_score(bench_only["y_true"], bench_only["CI_orig"])
        auc_final = roc_auc_score(bench_only["y_true"], bench_only["CI_final"])
        print(f"\nSubset AUC (Original): {auc_orig:.3f}")
        print(f"Subset AUC (Vetoed Final): {auc_final:.3f}")
        
    print("\nDetailed Comparison (Sorted by CI_final):")
    print(res_df.sort_values("CI_final", ascending=False)[["name", "label", "CI_orig", "SatRatio_pre", "CI_final"]])

if __name__ == "__main__":
    run_experiment()
