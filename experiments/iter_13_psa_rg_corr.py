
import sys
import os
import math
import numpy as np
from scipy.stats import pearsonr
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdFreeSASA

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))
import chameleon

def run_experiment():
    labels = {}
    with open("chameleon_local/labelled_set.tsv") as f:
        header = f.readline().strip().split("\t")
        iL = header.index("label")
        iN = header.index("name")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                labels[parts[iN]] = parts[iL]

    all_items = []
    with open("chameleon_local/benchmark.tsv") as f:
        for i, line in enumerate(f):
            if i > 12: break
            parts = line.strip().split("\t")
            if len(parts) == 2:
                all_items.append((parts[0], parts[1]))
    
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                all_items.append((parts[0], parts[1]))

    N_CONF = 15
    print(f"{'Name':<20} | {'Corr':>7} | {'Label'}")
    print("-" * 40)

    for name, smi in all_items:
        try:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xC0FFEE
            cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=N_CONF, params=params))
            AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500, mmffVariant="MMFF94s")
            
            radii = rdFreeSASA.classifyAtoms(mol)
            psas = []
            rgs = []
            for cid in cids:
                psas.append(chameleon.compute_3d_psa(mol, cid, radii))
                rgs.append(chameleon.compute_rg(mol, cid))
            
            if np.std(psas) < 1e-6 or np.std(rgs) < 1e-6:
                corr = 0.0
            else:
                corr, _ = pearsonr(psas, rgs)
            
            label = labels.get(name)
            is_pos = (label in ("chameleon", "protac")) or (name in ("protac_1", "protac_2"))
            l_str = "CHAM" if is_pos else "NON"
            
            print(f"{name[:20]:<20} | {corr:>7.3f} | {l_str}")
        except Exception as e:
            pass

if __name__ == "__main__":
    run_experiment()
