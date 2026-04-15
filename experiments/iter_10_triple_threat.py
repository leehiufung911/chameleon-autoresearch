
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))
import chameleon

def count_potential_imhb_pairs(mol):
    topo = Chem.GetDistanceMatrix(mol)
    donors = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in ("N", "O"):
            for n in a.GetNeighbors():
                if n.GetSymbol() == "H":
                    donors.append(a.GetIdx())
                    break
    acceptors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ("N", "O")]
    
    pairs = 0
    for d in donors:
        for a in acceptors:
            if d != a and topo[d, a] >= 4:
                pairs += 1
    return max(pairs, 1)

def get_max_alkyl_chain(mol):
    # Match consecutive non-ring CH2 groups
    for length in range(12, 1, -1):
        smarts = "[CH2X4;!R]" * length
        query = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(query):
            return length
    return 0

def run_experiment():
    user_mols = []
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                user_mols.append((parts[0], parts[1]))

    N_CONF = 20
    print(f"{'Name':<20} | {'IMHB':>5} | {'Nrot':>4} | {'Pot':>4} | {'Index':>7} | {'MaxAlk':>6}")
    print("-" * 65)

    for name, smi in user_mols:
        try:
            mol = Chem.MolFromSmiles(smi)
            mol_with_hs = Chem.AddHs(mol)
            summary = chameleon.summarize(name, smi, n_conf=N_CONF, verbose=False)
            
            imhb = summary.imhb_mean
            nrot = rdMolDescriptors.CalcNumRotatableBonds(mol)
            pot = count_potential_imhb_pairs(mol_with_hs)
            max_alk = get_max_alkyl_chain(mol)
            
            # Triple Threat Index
            # imhb * log(2+Nrot) / potential_sites
            index = imhb * math.log(2 + nrot) / pot
            
            print(f"{name:<20} | {imhb:>5.2f} | {nrot:>4} | {pot:>4} | {index:>7.4f} | {max_alk:>6}")
        except Exception as e:
            print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    run_experiment()
