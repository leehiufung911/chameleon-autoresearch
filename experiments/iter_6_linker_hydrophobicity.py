
import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Add the project root to sys.path to import from chameleon_local
sys.path.append(os.path.abspath("chameleon_local"))
import chameleon

def get_max_nonring_alkyl_chain(mol):
    # Match consecutive non-ring CH2 groups
    for length in range(10, 1, -1):
        smarts = "[CH2X4;!R]" * length
        query = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(query):
            return length
    return 0

def calculate_gated_alkyl_penalty(mol, imhb_mean):
    max_alk = get_max_nonring_alkyl_chain(mol)
    # Threshold: >= 4 carbons AND imhb_mean < 0.9
    if max_alk >= 4 and imhb_mean < 0.9:
        return 0.5
    return 1.0

def run_experiment():
    # Load user PROTACs
    user_mols = []
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                user_mols.append((parts[0], parts[1]))

    print(f"{'Name':<20} | {'CI_orig':>7} | {'MaxAlk':>6} | {'IMHB':>5} | {'CI_new':>7} | {'Verdict'}")
    print("-" * 75)

    N_CONF = 30 
    for name, smi in user_mols:
        try:
            summary = chameleon.summarize(name, smi, n_conf=N_CONF, verbose=False)
            mol = Chem.MolFromSmiles(smi)
            
            psa_range = summary.psa3d_max - summary.psa3d_min
            rg_ratio = summary.rg_max / max(summary.rg_min, 1e-6)
            term_psa = psa_range / math.sqrt(max(summary.mw, 1.0))
            term_imhb = summary.imhb_mean
            term_rg = math.log(max(rg_ratio, 1.0))
            
            # Coupled CI as base (from Iteration 6)
            ci_orig = term_psa + 10.0 * (term_imhb * term_rg)
            
            # Gated Alkyl Penalty
            penalty = calculate_gated_alkyl_penalty(mol, term_imhb)
            ci_new = ci_orig * penalty
            max_alk = get_max_nonring_alkyl_chain(mol)
            
            verdict = "CHAM" if ci_new > 4.0 else "NON" 
            print(f"{name:<20} | {ci_orig:>7.2f} | {max_alk:>6} | {term_imhb:>5.2f} | {ci_new:>7.2f} | {verdict}")
            
        except Exception as e:
            print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    run_experiment()
