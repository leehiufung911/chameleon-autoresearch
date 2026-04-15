
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.metrics import roc_auc_score

# Add chameleon_local to path
sys.path.append("chameleon_local")
import chameleon

def count_aliphatic_ethers(mol):
    # SMARTS for aliphatic ether: [O;H0;v2;!$(O=C)]
    # Matches Oxygen with no H, valence 2, not double bonded to Carbon
    query = Chem.MolFromSmarts("[O;H0;v2;!$(O=C)]")
    return len(mol.GetSubstructMatches(query))

def run_experiment():
    print("Loading data...", flush=True)
    # Benchmark
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    # Holdout
    df_hold = pd.read_csv("holdout_set.tsv", sep="\t")
    # User PROTACs
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    df_user["label"] = "user"

    # Sample for speed
    df_bench_sample = df_bench.sample(min(20, len(df_bench)), random_state=42)
    df_full = pd.concat([
        df_bench_sample[["name", "smiles", "label"]],
        df_hold[["name", "smiles", "label"]],
        df_user[["name", "smiles", "label"]]
    ], ignore_index=True)

    n_conf = 10
    results = []

    print(f"Starting processing of {len(df_full)} molecules...", flush=True)

    for i, row in df_full.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        if i % 10 == 0:
            print(f"  {i}/{len(df_full)}: {name}", flush=True)
        
        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None:
            continue
            
        mw = Descriptors.MolWt(mol0)
        ethers = count_aliphatic_ethers(mol0)
        
        # Conformer generation & MMFF minimization
        try:
            mol, cids = chameleon.embed_conformers(mol0, n_conf=n_conf, n_threads=0)
            energies = chameleon.minimize(mol, n_threads=0)
            # Prune to get the diversity
            keep, _ = chameleon.butina_prune(mol, energies)
        except Exception as e:
            print(f"Error processing {name}: {e}")
            continue

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
        
        # Raw descriptors
        psa_max, psa_min = psa.max(), psa.min()
        rg_max, rg_min = rg.max(), max(rg.min(), 0.1)
        imhb_raw = imhb.mean()
        
        # Components
        dPSA_norm = (psa_max - psa_min) / np.sqrt(mw)
        log_rg_ratio = np.log(rg_max / rg_min)
        
        # Formulas
        # 1. Iter 4 Winner: CI_mult
        ci_mult = (dPSA_norm + 3.0 * log_rg_ratio) * imhb_raw
        
        # 2. Pure Ratio: CI_ratio
        # We add 1.0 to IMHB to avoid 0s if we want to multiply, or just use additive
        ci_ratio = (psa_max / max(psa_min, 1.0)) * (rg_max / rg_min) * (imhb_raw + 0.1)
        
        # 3. Log Hybrid: (log(dPSA_ratio) + log(dRg_ratio)) * imhb
        # Equivalent to log(ratio_prod) * imhb
        ci_log_ratio = (np.log(psa_max / max(psa_min, 1.0)) + log_rg_ratio) * imhb_raw

        results.append({
            "name": name,
            "label": label,
            "mw": mw,
            "ethers": ethers,
            "imhb_raw": imhb_raw,
            "psa_ratio": psa_max / max(psa_min, 1.0),
            "rg_ratio": rg_max / rg_min,
            "ci_mult": ci_mult,
            "ci_ratio": ci_ratio,
            "ci_log_ratio": ci_log_ratio
        })

    df_res = pd.DataFrame(results)
    
    # Evaluation
    def eval_df(df_subset, title):
        print(f"\n--- {title} ---")
        df_eval = df_subset[df_subset["label"].isin(["chameleon", "nonchameleon"])].copy()
        if len(df_eval) < 2:
            print("Not enough data for AUC")
            return
        df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
        
        for col in ["imhb_raw", "ethers", "ci_mult", "ci_ratio", "ci_log_ratio"]:
            auc = roc_auc_score(df_eval["target"], df_eval[col])
            print(f"{col:15}: AUC = {auc:.3f}")

    # 1. Benchmark only
    eval_df(df_res[df_res["name"].isin(df_bench["name"])], "Benchmark Set (N=49)")
    
    # 2. Holdout only
    eval_df(df_res[df_res["name"].isin(df_hold["name"])], "Holdout Set (N=29)")
    
    # 3. Combined
    eval_df(df_res, "Combined Set (N=78)")

    print("\n--- User PROTACs ---")
    df_user_res = df_res[df_res["label"] == "user"].copy()
    print(df_user_res[["name", "ci_mult", "ci_ratio", "ci_log_ratio", "imhb_raw"]])

if __name__ == "__main__":
    try:
        run_experiment()
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)
