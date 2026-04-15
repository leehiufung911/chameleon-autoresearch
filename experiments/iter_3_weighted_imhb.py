
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
from sklearn.metrics import roc_auc_score

# Duplicate IMHB logic from chameleon.py but add ring size weighting
def get_weighted_imhb(mol, conf):
    donors = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ["N", "O"]:
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == "H":
                    donors.append((nbr.GetIdx(), atom.GetIdx()))
    acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ["N", "O"]]
    
    score = 0
    for h_idx, d_idx in donors:
        for a_idx in acceptors:
            if a_idx == d_idx: continue
            
            # Topological distance (D..A)
            path = Chem.GetShortestPath(mol, d_idx, a_idx)
            dist_topo = len(path) - 1 # number of bonds
            if dist_topo < 3: continue
            
            # Geometric check
            pos_h = np.array(conf.GetAtomPosition(h_idx))
            pos_a = np.array(conf.GetAtomPosition(a_idx))
            if np.linalg.norm(pos_h - pos_a) <= 2.5:
                pos_d = np.array(conf.GetAtomPosition(d_idx))
                v_hd = pos_d - pos_h
                v_ha = pos_a - pos_h
                angle = np.degrees(np.arccos(np.clip(np.dot(v_hd, v_ha) / (np.linalg.norm(v_hd) * np.linalg.norm(v_ha)), -1.0, 1.0)))
                if angle >= 120:
                    # WEIGHTING: 
                    # 5-10 bonds: 1.0 (ideal IMHB)
                    # > 10 bonds: 0.2 (accidental contact)
                    if 4 <= dist_topo <= 10:
                        score += 1.0
                    else:
                        score += 0.2
    return score

def run():
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    print("START", flush=True)
    results = []
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        print(f"Processing {name}...", flush=True)
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, numConfs=20, params=AllChem.ETKDGv3())
        # Minimize for realism
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
        
        imhb_weighted_total = 0
        imhb_raw_total = 0
        n = mol.GetNumConformers()
        for cid in range(n):
            conf = mol.GetConformer(cid)
            imhb_weighted_total += get_weighted_imhb(mol, conf)
            # Raw count
            imhb_raw_total += get_weighted_imhb(mol, conf) # Wait, I need a raw count function
            # Actually I'll just recompute it with weight 1.0
            
        # Re-get raw for comparison
        raw_score = 0
        for cid in range(n):
            conf = mol.GetConformer(cid)
            # Raw = weight 1.0 for all
            donors = [(nbr.GetIdx(), a.GetIdx()) for a in mol.GetAtoms() if a.GetSymbol() in ["N", "O"] for nbr in a.GetNeighbors() if nbr.GetSymbol()=="H"]
            acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ["N", "O"]]
            for h_idx, d_idx in donors:
                for a_idx in acceptors:
                    if a_idx == d_idx: continue
                    if len(Chem.GetShortestPath(mol, d_idx, a_idx)) <= 3: continue
                    pos_h, pos_a = np.array(conf.GetAtomPosition(h_idx)), np.array(conf.GetAtomPosition(a_idx))
                    if np.linalg.norm(pos_h - pos_a) <= 2.5:
                        pos_d = np.array(conf.GetAtomPosition(d_idx))
                        v_hd, v_ha = pos_d - pos_h, pos_a - pos_h
                        if np.degrees(np.arccos(np.clip(np.dot(v_hd, v_ha) / (np.linalg.norm(v_hd) * np.linalg.norm(v_ha)), -1.0, 1.0))) >= 120:
                            raw_score += 1
        
        results.append({
            "name": name,
            "imhb_raw": raw_score / n if n > 0 else 0,
            "imhb_weighted": imhb_weighted_total / n if n > 0 else 0
        })
        print(f"  {name}: Raw={results[-1]['imhb_raw']:.2f}, Weighted={results[-1]['imhb_weighted']:.2f}", flush=True)

    print("\n--- RESULTS ---")
    print(pd.DataFrame(results))

if __name__ == "__main__":
    run()
