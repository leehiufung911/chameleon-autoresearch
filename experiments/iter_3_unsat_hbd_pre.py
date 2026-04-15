
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
from sklearn.metrics import roc_auc_score

def get_imhb_pre(mol, n_conf=50):
    # Embed conformers
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=AllChem.ETKDGv3())
    if not cids:
        return 0
    
    # Geometric H-bond count (simplified version of chameleon.py)
    # H..A <= 2.5 A, D-H..A >= 120 deg
    
    # Get donors and acceptors
    # Using more precise definitions to match chameleon.py
    # D-H atoms
    donors = [] # list of (H_idx, D_idx)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ["N", "O"]:
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == "H":
                    donors.append((nbr.GetIdx(), atom.GetIdx()))
    
    # A atoms
    acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ["N", "O"]]
    
    imhb_total = 0
    for cid in cids:
        conf = mol.GetConformer(cid)
        count = 0
        for h_idx, d_idx in donors:
            for a_idx in acceptors:
                if a_idx == d_idx: continue
                # Topological distance > 3
                if len(Chem.GetShortestPath(mol, d_idx, a_idx)) <= 3: continue
                
                # Distance H..A
                dist = (np.array(conf.GetAtomPosition(h_idx)) - np.array(conf.GetAtomPosition(a_idx)))
                dist = np.linalg.norm(dist)
                if dist <= 2.5:
                    # Angle D-H..A
                    pos_d = np.array(conf.GetAtomPosition(d_idx))
                    pos_h = np.array(conf.GetAtomPosition(h_idx))
                    pos_a = np.array(conf.GetAtomPosition(a_idx))
                    v_dh = pos_h - pos_d
                    v_ha = pos_a - pos_h
                    # Cosine of angle
                    cos_theta = np.dot(v_dh, v_ha) / (np.linalg.norm(v_dh) * np.linalg.norm(v_ha))
                    # We want angle > 120 deg. Angle between v_dh and v_ha is 180 - theta.
                    # So theta should be small? No.
                    # v_dh and v_ha should be aligned. 
                    # Angle D-H-A is 180 - angle between (H-D) and (A-H).
                    # Actually, angle D-H..A is 180 - angle between v_hd and v_ha.
                    # Wait, let's use the standard formula.
                    # Angle is arccos(dot(v_hd, v_ha) / (|v_hd| * |v_ha|))
                    v_hd = pos_d - pos_h
                    v_ha = pos_a - pos_h
                    cos_val = np.dot(v_hd, v_ha) / (np.linalg.norm(v_hd) * np.linalg.norm(v_ha))
                    angle = np.degrees(np.arccos(np.clip(cos_val, -1.0, 1.0)))
                    if angle >= 120:
                        count += 1
        imhb_total += count
    
    return imhb_total / len(cids)

def run():
    print("Loading data...", flush=True)
    df_bench = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    N_CONF = 50
    results = []
    
    print(f"Processing {len(df_user)} User PROTACs...", flush=True)
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        total_hbd = rdMolDescriptors.CalcNumHBD(mol)
        imhb_pre = get_imhb_pre(mol, n_conf=N_CONF)
        unsat_pre = total_hbd - imhb_pre
        results.append({
            "name": name, "label": "protac", "hbd": total_hbd, "imhb_pre": imhb_pre, "unsat_pre": unsat_pre
        })
        print(f"  {name}: HBD={total_hbd}, IMHB_pre={imhb_pre:.2f}, Unsat_pre={unsat_pre:.2f}", flush=True)

    print(f"Processing Benchmark...", flush=True)
    for i, row in df_bench.iterrows():
        name, smiles, label = row["name"], row["smiles"], row["label"]
        mol = Chem.MolFromSmiles(smiles)
        total_hbd = rdMolDescriptors.CalcNumHBD(mol)
        imhb_pre = get_imhb_pre(mol, n_conf=N_CONF)
        unsat_pre = total_hbd - imhb_pre
        results.append({
            "name": name, "label": label, "hbd": total_hbd, "imhb_pre": imhb_pre, "unsat_pre": unsat_pre
        })
        if (i+1) % 10 == 0:
            print(f"  Progress: {i+1}/{len(df_bench)}", flush=True)

    df_res = pd.DataFrame(results)
    df_eval = df_res[df_res["label"].isin(["chameleon", "nonchameleon"])].copy()
    df_eval["target"] = (df_eval["label"] == "chameleon").astype(int)
    
    auc_imhb = roc_auc_score(df_eval["target"], df_eval["imhb_pre"])
    auc_unsat = roc_auc_score(df_eval["target"], -df_eval["unsat_pre"])
    
    print("\n--- AUC RESULTS (N_CONF=50, Pre-MMFF) ---", flush=True)
    print("AUC IMHB_pre:  {:.3f}".format(auc_imhb))
    print("AUC Unsat_pre: {:.3f}".format(auc_unsat))
    
    print("\n--- Summary Stats ---", flush=True)
    print(df_eval.groupby("label")[["imhb_pre", "unsat_pre"]].mean())
    
    print("\n--- User PROTACs Detail ---")
    print(df_res[df_res["label"] == "protac"][["name", "hbd", "imhb_pre", "unsat_pre"]])

if __name__ == "__main__":
    run()
