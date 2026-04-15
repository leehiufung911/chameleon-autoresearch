
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

def get_imhb_pre(mol, n_conf=10):
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=AllChem.ETKDGv3())
    if not cids: return 0
    donors = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ["N", "O"]:
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == "H":
                    donors.append((nbr.GetIdx(), atom.GetIdx()))
    acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ["N", "O"]]
    
    imhb_total = 0
    for cid in cids:
        conf = mol.GetConformer(cid)
        count = 0
        for h_idx, d_idx in donors:
            for a_idx in acceptors:
                if a_idx == d_idx: continue
                # Faster path check: just skip if topological distance <= 3
                if len(Chem.GetShortestPath(mol, d_idx, a_idx)) <= 3: continue
                pos_h = np.array(conf.GetAtomPosition(h_idx))
                pos_a = np.array(conf.GetAtomPosition(a_idx))
                if np.linalg.norm(pos_h - pos_a) <= 2.5:
                    pos_d = np.array(conf.GetAtomPosition(d_idx))
                    v_hd = pos_d - pos_h
                    v_ha = pos_a - pos_h
                    angle = np.degrees(np.arccos(np.clip(np.dot(v_hd, v_ha) / (np.linalg.norm(v_hd) * np.linalg.norm(v_ha)), -1.0, 1.0)))
                    if angle >= 120: count += 1
    return imhb_total / len(cids)

def run():
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    print("START", flush=True)
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        print(f"Processing {name}...", flush=True)
        mol = Chem.MolFromSmiles(smiles)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        imhb_pre = get_imhb_pre(mol, n_conf=20)
        print(f"  {name}: HBD={hbd}, IMHB_pre={imhb_pre:.4f}", flush=True)

if __name__ == "__main__":
    run()
