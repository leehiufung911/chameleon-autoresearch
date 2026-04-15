
import sys
import os

# Fix for Windows environment issues when running via direct python.exe
dll_path = r'C:\Users\mic23\miniconda3\envs\chameleon\Library\bin'
if os.path.exists(dll_path):
    os.environ['PATH'] = dll_path + os.pathsep + os.environ['PATH']
    if hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(dll_path)

print("Imports start...", flush=True)
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.metrics import roc_auc_score

sys.path.append("chameleon_local")
import chameleon
print("Imports done.", flush=True)

def get_hbd_count(mol):
    return rdMolDescriptors.CalcNumHBD(mol)

def compute_custom_imhb(mol, conf_id, topo_dist, min_topo=4):
    conf = mol.GetConformer(conf_id)
    coords = np.asarray(conf.GetPositions())
    donor_pairs = []
    acceptors = []
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        if sym in ("N", "O"):
            acceptors.append(a.GetIdx())
            for n in a.GetNeighbors():
                if n.GetSymbol() == "H":
                    donor_pairs.append((a.GetIdx(), n.GetIdx()))
    if not donor_pairs or not acceptors:
        return 0, 0
    count_all = 0
    count_global = 0
    used_donors = set()
    d_max = 2.5
    ang_min_rad = np.radians(120.0)
    for d_idx, h_idx in donor_pairs:
        if h_idx in used_donors: continue
        for a_idx in acceptors:
            if a_idx == d_idx: continue
            dist_2d = topo_dist[d_idx, a_idx]
            if dist_2d < min_topo: continue
            d_ha = np.linalg.norm(coords[h_idx] - coords[a_idx])
            if d_ha > d_max: continue
            v1 = coords[d_idx] - coords[h_idx]
            v2 = coords[a_idx] - coords[h_idx]
            nv1 = np.linalg.norm(v1)
            nv2 = np.linalg.norm(v2)
            if nv1 < 1e-6 or nv2 < 1e-6: continue
            cosang = float(np.dot(v1, v2) / (nv1 * nv2))
            cosang = max(-1.0, min(1.0, cosang))
            ang = np.arccos(cosang)
            if ang >= ang_min_rad:
                count_all += 1
                if dist_2d > 8: count_global += 1
                used_donors.add(h_idx)
                break
    return count_all, count_global

def run_experiment():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    df_bench_sample = df_bench.sample(min(20, len(df_bench)), random_state=42)
    # Filter out giants for speed
    df_bench_sample = df_bench_sample[df_bench_sample['smiles'].apply(lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x)) < 1200)]
    
    df_user["label"] = "user"
    df_full = pd.concat([df_bench_sample, df_user], ignore_index=True)
    
    results = []
    n_conf = 20
    n_threads = 1 # Use single thread to avoid hangs
    
    for i, row in df_full.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        print(f"Processing {name} ({i+1}/{len(df_full)})...", flush=True)
        
        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None: 
            print(f"  Failed to parse SMILES for {name}", flush=True)
            continue
        mw = Descriptors.MolWt(mol0)
        hbd = max(1, get_hbd_count(mol0))
        
        # 1. ETKDG (Pre-minimization) for SatRatio
        print(f"  Embedding ETKDG...", flush=True)
        mol_pre, cids_pre = chameleon.embed_conformers(mol0, n_conf=n_conf, n_threads=n_threads)
        topo = chameleon.bond_path_distance_matrix(mol_pre)
        
        imhb_pre_list = []
        for cid in cids_pre:
            c_all, _ = compute_custom_imhb(mol_pre, cid, topo)
            imhb_pre_list.append(c_all)
        
        sat_ratio_pre = np.mean(imhb_pre_list) / hbd
        
        # 2. MMFF (Post-minimization)
        print(f"  Embedding MMFF...", flush=True)
        mol_mmff, cids_mmff = chameleon.embed_conformers(mol0, n_conf=n_conf, n_threads=n_threads)
        print(f"  Minimizing...", flush=True)
        energies = chameleon.minimize(mol_mmff, n_threads=n_threads)
        print(f"  Pruning...", flush=True)
        keep, _ = chameleon.butina_prune(mol_mmff, energies)
        
        psa_list, rg_list, imhb_list, imhb_glob_list = [], [], [], []
        radii = chameleon.rdFreeSASA.classifyAtoms(mol_mmff)
        topo = chameleon.bond_path_distance_matrix(mol_mmff)
        
        for cid in keep:
            psa_list.append(chameleon.compute_3d_psa(mol_mmff, cid, radii))
            rg_list.append(chameleon.compute_rg(mol_mmff, cid))
            c_all, c_glob = compute_custom_imhb(mol_mmff, cid, topo)
            imhb_list.append(c_all)
            imhb_glob_list.append(c_glob)
            
        psa = np.array(psa_list)
        rg = np.array(rg_list)
        imhb = np.array(imhb_list)
        imhb_glob = np.array(imhb_glob_list)
        
        psa_range = psa.max() - psa.min()
        rg_ratio = rg.max() / max(rg.min(), 1e-6)
        
        term_psa = psa_range / np.sqrt(mw)
        term_rg = np.log(max(rg_ratio, 1.0))
        term_imhb = np.mean(imhb)
        term_imhb_glob = np.mean(imhb_glob)
        
        ci_orig = 2.0 * term_psa + 1.0 * term_imhb + 3.0 * term_rg
        ci_global = 2.0 * term_psa + 1.0 * term_imhb_glob + 3.0 * term_rg
        ci_coupled = term_psa + 10.0 * (term_imhb * term_rg)
        veto_005 = 1.0 if sat_ratio_pre > 0.005 else 0.4
        
        results.append({
            "name": name, "label": label, "ci_orig": ci_orig, "ci_global": ci_global,
            "ci_coupled": ci_coupled, "ci_vetoed_005": ci_orig * veto_005,
            "sat_ratio_pre": sat_ratio_pre, "imhb_mean": term_imhb,
            "imhb_glob_mean": term_imhb_glob, "term_psa": term_psa, "term_rg": term_rg
        })

    df_res = pd.DataFrame(results)
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    if not df_eval.empty:
        df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
        print("\n--- AUC Results (Subsample) ---", flush=True)
        for col in ["ci_orig", "ci_global", "ci_coupled", "ci_vetoed_005"]:
            print(f"{col:15}: {roc_auc_score(df_eval['target'], df_eval[col]):.3f}", flush=True)
        
    print("\n--- User PROTACs ---", flush=True)
    print(df_res[df_res["label"] == "user"][["name", "ci_orig", "ci_vetoed_005", "ci_coupled", "sat_ratio_pre"]], flush=True)

if __name__ == "__main__":
    run_experiment()
