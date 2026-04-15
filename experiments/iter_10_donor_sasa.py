
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFreeSASA
import csv

# Add chameleon_local to path
sys.path.append("chameleon_local")
import chameleon

def get_donor_sasa(mol, conf_id, radii):
    rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id)
    donor_indices = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in ("N", "O"):
            if any(n.GetSymbol() == "H" for n in a.GetNeighbors()):
                donor_indices.append(a.GetIdx())
                # and the H atoms themselves
                for n in a.GetNeighbors():
                    if n.GetSymbol() == "H":
                        donor_indices.append(n.GetIdx())
    
    if not donor_indices:
        return 0.0
        
    sasa = 0.0
    for idx in donor_indices:
        sasa += float(mol.GetAtomWithIdx(idx).GetProp("SASA"))
    return sasa

def run_experiment():
    items = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) == 2: items.append((row[0], row[1]))

    bench_subset = ["cyclosporinA", "tacrolimus", "aspirin", "ibuprofen"]
    with open("chameleon_local/benchmark.tsv", "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) == 2 and row[0] in bench_subset:
                items.append((row[0], row[1]))

    N_CONF = 30
    print(f"{'Name':<20} | {'D_min':>7} | {'D_max':>7} | {'D_ratio':>7}")
    print("-" * 55)

    for name, smi in items:
        try:
            mol0 = Chem.MolFromSmiles(smi)
            mol, cids = chameleon.embed_conformers(mol0, n_conf=N_CONF, n_threads=1)
            chameleon.minimize(mol, n_threads=1)
            radii = rdFreeSASA.classifyAtoms(mol)
            
            sasa_list = [get_donor_sasa(mol, cid, radii) for cid in range(mol.GetNumConformers())]
            sasa_list = [s for s in sasa_list if s > 0]
            if not sasa_list:
                print(f"{name:<20} | NO DONORS")
                continue
                
            d_min = min(sasa_list)
            d_max = max(sasa_list)
            d_ratio = d_min / d_max
            
            print(f"{name:<20} | {d_min:>7.1f} | {d_max:>7.1f} | {d_ratio:>7.3f}")
            
        except Exception as e:
            print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    run_experiment()
