
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))
import chameleon

def get_sa_metrics(mol, confId):
    # RDKit's LabuteASA is a 2D descriptor, doesn't take confId.
    # We need a 3D surface area. 
    # Let's use rdFreeSASA if available, or just a simple VdW surface approximation.
    # Actually, chameleon.compute_psa3d uses an approximation.
    # Let's use a similar approach for total surface area.
    
    psa3d = chameleon.compute_psa3d(mol, confId)
    
    # Simple 3D SA approximation: sum of VdW areas of atoms
    # Or just use the fact that we have 3D coords.
    # Let's use a simple sphere-sum with overlap correction? Too complex.
    # What if we use rdMolDescriptors.CalcTPSA for PSA and something else for Total?
    # Actually, chameleon's psa3d is the best we have for 3D PSA.
    # For total 3D SA, we can use a similar heuristic.
    
    conf = mol.GetConformer(confId)
    total_sa = 0
    # VdW radii (approx)
    radii = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "F": 1.47, "Cl": 1.75, "Br": 1.85, "P": 1.8, "H": 1.2}
    
    for i in range(mol.GetNumAtoms()):
        symbol = mol.GetAtomWithIdx(i).GetSymbol()
        r = radii.get(symbol, 1.5)
        total_sa += 4 * math.pi * (r**2)
    # This is constant across conformers if we don't account for overlap!
    # We MUST account for overlap to see the difference between compact and extended.
    
    # Let's use the Rg as a proxy for total surface area?
    # SA ~ Rg^2. 
    rg = chameleon.compute_rg(mol, confId)
    return psa3d, rg

def run_experiment():
    user_mols = []
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                user_mols.append((parts[0], parts[1]))

    N_CONF = 20
    print(f"{'Name':<20} | {'dPSA':>6} | {'dRg':>6} | {'Ratio':>7}")
    print("-" * 50)

    for name, smi in user_mols:
        try:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xC0FFEE
            cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=N_CONF, params=params))
            AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500, mmffVariant="MMFF94s")
            
            psas = []
            rgs = []
            for cid in cids:
                psas.append(chameleon.compute_psa3d(mol, cid))
                rgs.append(chameleon.compute_rg(mol, cid))
            
            d_psa = np.max(psas) - np.min(psas)
            d_rg = np.max(rgs) - np.min(rgs)
            
            # Ratio of folding to masking
            # If a molecule collapses a lot (dRg) but doesn't mask much PSA (dPSA), 
            # then Ratio (dRg/dPSA) is high. 
            # This indicates "inefficient" or "hydrophobic" collapse.
            ratio = (d_rg / d_psa * 100) if d_psa > 0 else 0
            
            print(f"{name:<20} | {d_psa:>6.1f} | {d_rg:>6.2f} | {ratio:>7.3f}")
        except Exception as e:
            print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    run_experiment()
