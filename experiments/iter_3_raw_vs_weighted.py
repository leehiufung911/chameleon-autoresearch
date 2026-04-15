
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import roc_auc_score

sys.path.append("chameleon_local")
import chameleon

def get_imhb(mol, conf):
    donors = [(nbr.GetIdx(), a.GetIdx()) for a in mol.GetAtoms() if a.GetSymbol() in ["N", "O"] for nbr in a.GetNeighbors() if nbr.GetSymbol()=="H"]
    acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ["N", "O"]]
    count = 0
    for h_idx, d_idx in donors:
        for a_idx in acceptors:
            if a_idx == d_idx: continue
            if len(Chem.GetShortestPath(mol, d_idx, a_idx)) <= 3: continue
            pos_h, pos_a = np.array(conf.GetAtomPosition(h_idx)), np.array(conf.GetAtomPosition(a_idx))
            if np.linalg.norm(pos_h - pos_a) <= 2.5:
                pos_d = np.array(conf.GetAtomPosition(d_idx))
                v_hd, v_ha = pos_d - pos_h, pos_a - pos_h
                if np.degrees(np.arccos(np.clip(np.dot(v_hd, v_ha) / (np.linalg.norm(v_hd) * np.linalg.norm(v_ha)), -1.0, 1.0))) >= 120:
                    count += 1
    return count

def run():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    # Sample for speed (10 each)
    df_eval_set = pd.concat([df_bench[df_bench["label"]=="chameleon"].sample(10, random_state=42), 
                             df_bench[df_bench["label"]=="nonchameleon"].sample(10, random_state=42)])
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    N_CONF = 20
    results = []
    
    for df, label_type in [(df_user, "protac"), (df_eval_set, "benchmark")]:
        print(f"Processing {label_type}...", flush=True)
        for i, row in df.iterrows():
            name, smiles = row["name"], row["smiles"]
            label = row["label"] if "label" in row else "protac"
            try:
                mol = Chem.MolFromSmiles(smiles)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMultipleConfs(mol, numConfs=N_CONF, params=AllChem.ETKDGv3())
                energies = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
                
                imhbs = []
                n = mol.GetNumConformers()
                for cid in range(n):
                    imhbs.append(get_imhb(mol, mol.GetConformer(cid)))
                
                # Raw mean
                imhb_raw = np.mean(imhbs)
                
                # Weighted mean (Boltzmann proxy like chameleon.py)
                # E_eff = E_mmff - 2.0 * n_imhb
                e_mmff = np.array([e[1] for e in energies])
                e_eff = e_mmff - 2.0 * np.array(imhbs)
                e_rel = e_eff - np.min(e_eff)
                weights = np.exp(-e_rel / 0.5924)
                weights /= np.sum(weights)
                imhb_weighted = np.sum(np.array(imhbs) * weights)
                
                results.append({
                    "name": name, "label": label, "raw": imhb_raw, "weighted": imhb_weighted
                })
                if label_type == "protac":
                    print(f"  {name}: Raw={imhb_raw:.2f}, Weighted={imhb_weighted:.2f}", flush=True)
            except: pass

    df_res = pd.DataFrame(results)
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
    
    print("\n--- AUC (n=20) ---")
    print("AUC Raw:      {:.3f}".format(roc_auc_score(df_eval["target"], df_eval["raw"])))
    print("AUC Weighted: {:.3f}".format(roc_auc_score(df_eval["target"], df_eval["weighted"])))
    
    print("\n--- User PROTACs Detail ---")
    print(df_res[df_res["label"] == "protac"][["name", "raw", "weighted"]])

if __name__ == "__main__":
    run()
