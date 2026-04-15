import sys
import os
import numpy as np
from scipy.stats import pearsonr
from rdkit import Chem

sys.path.append(os.path.abspath("chameleon_local"))
from chameleon import (
    embed_conformers, minimize, butina_prune, 
    compute_3d_psa, compute_rg, rdFreeSASA
)

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
                    
    print(f"{'Name':<15} | {'Label':<5} | {'Pearson(Rg,PSA)':>15} | {'p-value':>8}")
    print("-" * 55)
    
    for name, smi in all_items:
        mol0 = Chem.MolFromSmiles(smi)
        if not mol0: continue
        
        mol, cids = embed_conformers(mol0, 100, 0)
        energies = minimize(mol, 0)
        keep, _ = butina_prune(mol, energies, rms_cut=0.75)
        
        radii = rdFreeSASA.classifyAtoms(mol)
        psa_list, rg_list = [], []
        for cid in keep:
            psa_list.append(compute_3d_psa(mol, cid, radii))
            rg_list.append(compute_rg(mol, cid))
            
        if len(psa_list) > 2:
            corr, pval = pearsonr(rg_list, psa_list)
            is_pos = targets[name]
            l_str = "CHAM" if is_pos else "NON"
            print(f"{name:<15} | {l_str:<5} | {corr:15.3f} | {pval:8.3f}")

if __name__ == "__main__":
    try:
        run_experiment()
    except Exception as e:
        import traceback
        traceback.print_exc()
