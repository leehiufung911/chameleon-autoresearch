
import os
import sys

# Fix for Windows environment issues when running via direct python.exe
dll_path = r'C:\Users\mic23\miniconda3\envs\chameleon\Library\bin'
if os.path.exists(dll_path):
    os.environ['PATH'] = dll_path + os.pathsep + os.environ['PATH']
    if hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(dll_path)

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors
import math

def get_imhb_count(mol, conf_id, topo_dist, d_max=2.5, ang_min_deg=120.0, min_topo=4):
    conf = mol.GetConformer(conf_id)
    coords = np.asarray(conf.GetPositions())
    ang_min_rad = math.radians(ang_min_deg)
    donor_pairs = []
    acceptors = []
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        if sym in ("N", "O"):
            acceptors.append(a.GetIdx())
            for n in a.GetNeighbors():
                if n.GetSymbol() == "H":
                    donor_pairs.append((a.GetIdx(), n.GetIdx()))
    if not donor_pairs or not acceptors: return 0
    count = 0
    used_donors = set()
    for d_idx, h_idx in donor_pairs:
        if h_idx in used_donors: continue
        for a_idx in acceptors:
            if a_idx == d_idx: continue
            if topo_dist[d_idx, a_idx] < min_topo: continue
            d_ha = np.linalg.norm(coords[h_idx] - coords[a_idx])
            if d_ha > d_max: continue
            v1, v2 = coords[d_idx] - coords[h_idx], coords[a_idx] - coords[h_idx]
            nv1, nv2 = np.linalg.norm(v1), np.linalg.norm(v2)
            if nv1 < 1e-6 or nv2 < 1e-6: continue
            cosang = max(-1.0, min(1.0, np.dot(v1, v2) / (nv1 * nv2)))
            if math.acos(cosang) >= ang_min_rad:
                count += 1
                used_donors.add(h_idx)
                break
    return count

def compute_3d_psa(mol, conf_id):
    from rdkit.Chem import rdFreeSASA
    radii = rdFreeSASA.classifyAtoms(mol)
    rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id)
    
    psa = 0.0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ('N', 'O'):
            psa += float(atom.GetProp("SASA"))
    return psa

def run_experiment():
    benchmark = pd.read_csv('chameleon_local/labelled_set.tsv', sep='\t')
    user_protacs = pd.read_csv('chameleon_local/user_protacs.tsv', sep='\t', names=['name', 'smiles'])
    user_protacs['label'] = [1, 1, 0]
    benchmark['label'] = benchmark['label'].apply(lambda x: 1 if x in ['chameleon', 'protac'] else 0)
    combined = pd.concat([benchmark, user_protacs], ignore_index=True)
    
    results = []
    for idx, row in combined.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol is None: continue
        if mol.GetNumHeavyAtoms() > 120: continue # Faster
        
        print(f"Processing {row['name']}...", flush=True)
        mol = Chem.AddHs(mol)
        topo_dist = Chem.GetDistanceMatrix(mol)
        AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=42) # 10 confs for speed
        AllChem.MMFFOptimizeMoleculeConfs(mol)
        
        imhb_counts = []
        psa_3d_values = []
        for i in range(mol.GetNumConformers()):
            imhb_counts.append(get_imhb_count(mol, i, topo_dist))
            psa_3d_values.append(compute_3d_psa(mol, i))
            
        imhb_raw = np.mean(imhb_counts)
        psa_ratio = np.max(psa_3d_values) / np.min(psa_3d_values) if np.min(psa_3d_values) > 0 else 1.0
        psa_diff = np.max(psa_3d_values) - np.min(psa_3d_values)
        mw = Descriptors.MolWt(mol)
        n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        
        results.append({
            'name': row['name'],
            'label': row['label'],
            'imhb_raw': imhb_raw,
            'psa_ratio': psa_ratio,
            'psa_diff_norm': psa_diff / math.sqrt(mw),
            'imhb_density': imhb_raw / (mw / 100.0),
            'imhb_rot': imhb_raw * math.log1p(n_rot)
        })
        
    df = pd.DataFrame(results)
    from sklearn.metrics import roc_auc_score
    
    print("\nRESULTS (AUC):")
    for col in ['imhb_raw', 'psa_ratio', 'psa_diff_norm', 'imhb_density', 'imhb_rot']:
        auc = roc_auc_score(df['label'], df[col])
        print(f"AUC {col:<15}: {auc:.3f}")
        
    # Hybrid metrics
    df['hybrid_ratio'] = np.log(df['psa_ratio']) * df['imhb_raw']
    df['hybrid_diff'] = df['psa_diff_norm'] * df['imhb_raw']
    
    print(f"AUC hybrid_ratio   : {roc_auc_score(df['label'], df['hybrid_ratio']):.3f}")
    print(f"AUC hybrid_diff    : {roc_auc_score(df['label'], df['hybrid_diff']):.3f}")

    print("\nPROTACs:")
    protacs = df[df['name'].str.contains('protac')].copy()
    print(protacs[['name', 'label', 'imhb_raw', 'psa_ratio', 'hybrid_ratio']])

if __name__ == "__main__":
    run_experiment()
