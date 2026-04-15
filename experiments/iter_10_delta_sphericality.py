import sys
import os
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

sys.path.append(os.path.abspath("chameleon_local"))
from chameleon import (
    embed_conformers, minimize, butina_prune, 
    compute_3d_psa, compute_rg, rdFreeSASA,
    compute_imhb, bond_path_distance_matrix,
    summarize
)

def analyze_shape(name, smiles, n_conf=100):
    mol0 = Chem.MolFromSmiles(smiles)
    if not mol0:
        return None
    
    # 1. Generate & prune conformers
    mol, cids = embed_conformers(mol0, n_conf, 0)
    energies = minimize(mol, 0)
    keep, k_energies = butina_prune(mol, energies, rms_cut=0.75)
    
    # 2. Compute descriptors
    radii = rdFreeSASA.classifyAtoms(mol)
    
    psa_list, rg_list, asph_list = [], [], []
    for cid in keep:
        psa_list.append(compute_3d_psa(mol, cid, radii))
        rg_list.append(compute_rg(mol, cid))
        asph_list.append(rdMolDescriptors.CalcAsphericity(mol, confId=cid))
        
    psa = np.asarray(psa_list)
    rg = np.asarray(rg_list)
    asph = np.asarray(asph_list)
    
    if len(psa) == 0:
        return None
        
    # 3. Define "polar" (extended, high PSA) vs "apolar" (compact, low PSA) sub-ensembles
    # Using the bottom/top 20% by PSA as in the main script
    q_size = max(1, int(len(psa) * 0.2))
    
    idx_sorted = np.argsort(psa)
    idx_compact = idx_sorted[:q_size]
    idx_extended = idx_sorted[-q_size:]
    
    asph_compact = np.mean(asph[idx_compact])
    asph_extended = np.mean(asph[idx_extended])
    delta_asph = asph_extended - asph_compact
    
    return {
        "asph_min": np.min(asph),
        "asph_max": np.max(asph),
        "asph_compact": asph_compact,
        "asph_extended": asph_extended,
        "delta_asph": delta_asph
    }

def run_experiment():
    targets = {
        "protac_1": True,
        "protac_2": True,
        "protac_3": False,
        "dBET6": True,
        "cyclosporinA": True,
        "atorvastatin": False
    }
    
    all_items = []
    for tsv in ["chameleon_local/benchmark.tsv", "chameleon_local/user_protacs.tsv"]:
        with open(tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts[0] in targets:
                    all_items.append((parts[0], parts[1]))
                    
    print(f"{'Name':<15} | {'Label':<5} | {'Asph(Comp)':>10} | {'Asph(Ext)':>9} | {'DeltaAsph':>9}")
    print("-" * 60)
    
    for name, smi in all_items:
        res = analyze_shape(name, smi, n_conf=100)
        if res:
            is_pos = targets[name]
            l_str = "CHAM" if is_pos else "NON"
            
            # Asphericity: 0 = perfect sphere, 1 = perfect rod
            # Delta = Extended - Compact. 
            # We expect true chameleons to be MORE spherical (lower asph) when compact.
            # So DeltaAsph should be positive and large for chameleons.
            
            print(f"{name:<15} | {l_str:<5} | {res['asph_compact']:10.3f} | {res['asph_extended']:9.3f} | {res['delta_asph']:9.3f}")

if __name__ == "__main__":
    try:
        run_experiment()
    except Exception as e:
        import traceback
        traceback.print_exc()
