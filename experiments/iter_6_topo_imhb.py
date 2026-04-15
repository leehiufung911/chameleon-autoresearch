
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
from rdkit.Chem import AllChem, rdMolDescriptors
import math
import sys

def get_imhb_data(mol, conf_id, topo_dist, d_max=2.5, ang_min_deg=120.0, min_topo=4):
    """
    Returns a list of topological distances for all IMHBs found in the conformer.
    Uses the same geometric criteria as chameleon.py
    """
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
    
    if not donor_pairs or not acceptors:
        return []

    imhb_distances = []
    used_donors = set()
    for d_idx, h_idx in donor_pairs:
        if h_idx in used_donors:
            continue
        for a_idx in acceptors:
            if a_idx == d_idx:
                continue
            if topo_dist[d_idx, a_idx] < min_topo:
                continue
            
            # Distance H..A
            d_ha = np.linalg.norm(coords[h_idx] - coords[a_idx])
            if d_ha > d_max:
                continue
            
            # Angle D-H..A
            v1 = coords[d_idx] - coords[h_idx]
            v2 = coords[a_idx] - coords[h_idx]
            nv1 = np.linalg.norm(v1)
            nv2 = np.linalg.norm(v2)
            if nv1 < 1e-6 or nv2 < 1e-6:
                continue
            cosang = float(np.dot(v1, v2) / (nv1 * nv2))
            cosang = max(-1.0, min(1.0, cosang))
            ang = math.acos(cosang)
            
            if ang >= ang_min_rad:
                imhb_distances.append(int(topo_dist[d_idx, a_idx]))
                used_donors.add(h_idx)
                break
                
    return imhb_distances

def run_experiment():
    # Load benchmark and labels
    benchmark = pd.read_csv('chameleon_local/labelled_set.tsv', sep='\t')
    
    # Add user PROTACs
    user_protacs = pd.read_csv('chameleon_local/user_protacs.tsv', sep='\t', names=['name', 'smiles'])
    user_protacs['label'] = [1, 1, 0] # 1, 2 are chameleonic, 3 is not
    
    # Map labels: chameleon and protac are positive (1), nonchameleon is negative (0)
    benchmark['label'] = benchmark['label'].apply(lambda x: 1 if x in ['chameleon', 'protac'] else 0)
    
    combined = pd.concat([benchmark, user_protacs], ignore_index=True)
    
    results = []
    
    print(f"Processing {len(combined)} molecules...")
    
    for idx, row in combined.iterrows():
        smiles = row['smiles']
        name = row['name']
        label = row['label']
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: continue
        
        # Skip very large molecules to stay under 5 mins
        if mol.GetNumHeavyAtoms() > 150:
            print(f"Skipping {name} (too large: {mol.GetNumHeavyAtoms()} atoms)")
            continue
            
        print(f"Processing {name} ({mol.GetNumHeavyAtoms()} heavy atoms)...", flush=True)
        mol = Chem.AddHs(mol)
        
        topo_dist = Chem.GetDistanceMatrix(mol)
        
        # Generate ensemble
        res = AllChem.EmbedMultipleConfs(mol, numConfs=20, randomSeed=42, pruneRmsThresh=0.5)
        if not res:
            # Fallback if embedding fails
            AllChem.EmbedMultipleConfs(mol, numConfs=20, randomSeed=42, useRandomCoords=True)
            
        AllChem.MMFFOptimizeMoleculeConfs(mol)
        
        all_imhb_dists = []
        n_confs = mol.GetNumConformers()
        for i in range(n_confs):
            dists = get_imhb_data(mol, i, topo_dist)
            all_imhb_dists.append(dists)
            
        # Descriptors
        # 1. Raw IMHB count
        counts = [len(d) for d in all_imhb_dists]
        imhb_raw = np.mean(counts) if counts else 0
        
        # 2. Long-range IMHB count (> 8 bonds)
        long_counts = [len([dist for dist in d if dist > 8]) for d in all_imhb_dists]
        imhb_long = np.mean(long_counts) if long_counts else 0
        
        # 3. Entropy-penalized IMHB (sum 1/log(dist))
        # Note: log(4) is the min distance. We use 1/log2(dist) maybe?
        # Let's use 1/log(dist)
        def entropy_score(dists):
            return sum(1.0 / math.log(d) for d in dists)
        
        entropy_scores = [entropy_score(d) for d in all_imhb_dists]
        imhb_entropy = np.mean(entropy_scores) if entropy_scores else 0
        
        results.append({
            'name': name,
            'label': label,
            'imhb_raw': imhb_raw,
            'imhb_long': imhb_long,
            'imhb_entropy': imhb_entropy
        })
        
    df = pd.DataFrame(results)
    
    from sklearn.metrics import roc_auc_score
    
    auc_raw = roc_auc_score(df['label'], df['imhb_raw'])
    auc_long = roc_auc_score(df['label'], df['imhb_long'])
    auc_entropy = roc_auc_score(df['label'], df['imhb_entropy'])
    
    print(f"\nRESULTS (N=20 conformers):")
    print(f"AUC IMHB Raw:     {auc_raw:.3f}")
    print(f"AUC IMHB Long:    {auc_long:.3f}")
    print(f"AUC IMHB Entropy: {auc_entropy:.3f}")
    
    print("\nPROTACs:")
    protacs = df[df['name'].str.contains('protac')].copy()
    print(protacs[['name', 'label', 'imhb_raw', 'imhb_long', 'imhb_entropy']])

if __name__ == "__main__":
    print("Starting script...")
    run_experiment()
