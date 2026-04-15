import os
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
import pandas as pd
import time

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
try:
    from chameleon_local.chameleon import compute_3d_psa, get_radii
except ImportError:
    # If chameleon.py cannot be imported, we define a dummy or duplicate
    pass

def compute_psa(mol, conf_id, radii):
    try:
        from chameleon_local.chameleon import compute_3d_psa
        return compute_3d_psa(mol, conf_id, radii)
    except:
        return rdMolDescriptors.CalcTPSA(mol)

def run_experiment():
    # Small test set: 3 chameleons, 1 non-chameleon, plus user PROTACs
    smiles_dict = {
        "cyclosporinA": "CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1C)C(C)C)C)CC2=CC=CC=C2)C)C)C)CC(C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C(C)C", # Just simplified
        "protac_1": "O=C(C(N1C(C2=CC=CC(NCCOCCOCCOCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=C2C1=O)=O)CC7)NC7=O",
        "protac_2": "O=C(C(N1C(C2=CC=CC(NC(COCCOCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
        "protac_3": "O=C(C(N1C(C2=CC=CC(OCC(NCCCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
        "atorvastatin": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC(O)CC(O)CC(=O)O",
    }
    
    # We will load actual SMILES from files to be safe
    molecules = []
    
    # Read user PROTACs
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            if not line.strip(): continue
            parts = line.strip().split()
            if len(parts) >= 2:
                molecules.append((parts[0], parts[1], "CHAM" if "1" in parts[0] or "2" in parts[0] else "NON"))
                
    # Add a few benchmarks
    bench_found = 0
    with open("chameleon_local/labelled_set.tsv") as f:
        for line in f:
            if bench_found >= 3: break
            if "cyclosporin" in line.lower() or "atorvastatin" in line.lower() or "tacrolimus" in line.lower():
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    molecules.append((parts[0], parts[1], parts[2]))
                    bench_found += 1
                    
    print(f"{'Name':<15} {'Label':<5} {'dPSA':<8} {'dAsph':<8} {'dSph':<8} {'Asph_min':<8} {'Asph_max':<8}")
    
    for name, smiles, label in molecules:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: continue
        mol = Chem.AddHs(mol)
        
        n_conf = 50
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        
        AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
        
        if mol.GetNumConformers() == 0:
            continue
            
        # MMFF Minimize
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
        
        try:
            from chameleon_local.chameleon import get_radii
            radii = get_radii()
        except:
            radii = None
            
        data = []
        for conf in mol.GetConformers():
            cid = conf.GetId()
            psa = compute_psa(mol, cid, radii)
            asph = rdMolDescriptors.CalcAsphericity(mol, confId=cid)
            sph = rdMolDescriptors.CalcSpherocityIndex(mol, confId=cid)
            data.append((psa, asph, sph))
            
        data.sort(key=lambda x: x[0]) # sort by PSA
        
        q = max(1, len(data) // 5)
        apolar = data[:q]
        polar = data[-q:]
        
        psa_min = np.mean([x[0] for x in apolar])
        psa_max = np.mean([x[0] for x in polar])
        asph_min_psa = np.mean([x[1] for x in apolar])
        asph_max_psa = np.mean([x[1] for x in polar])
        sph_min_psa = np.mean([x[2] for x in apolar])
        sph_max_psa = np.mean([x[2] for x in polar])
        
        dPSA = psa_max - psa_min
        dAsph = asph_max_psa - asph_min_psa
        dSph = sph_max_psa - sph_min_psa
        
        print(f"{name:<15} {label:<5} {dPSA:<8.1f} {dAsph:<8.3f} {dSph:<8.3f} {asph_min_psa:<8.3f} {asph_max_psa:<8.3f}")

if __name__ == "__main__":
    t0 = time.time()
    run_experiment()
    print(f"Time: {time.time()-t0:.1f}s")
