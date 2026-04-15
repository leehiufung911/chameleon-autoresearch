import sys
import os
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

sys.path.append(os.path.abspath("chameleon_local"))
import chameleon

def get_shape_metrics(mol, keep, w_apolar, w_polar):
    asph_list = []
    for cid in keep:
        asph = rdMolDescriptors.CalcAsphericity(mol, confId=cid)
        asph_list.append(asph)
    
    asph = np.asarray(asph_list)
    asph_ap = float(np.sum(w_apolar * asph))
    asph_po = float(np.sum(w_polar * asph))
    
    return asph_ap, asph_po, asph.min(), asph.max(), asph.mean()

def run_experiment():
    targets = {
        "cyclosporinA": True,
        "tacrolimus": True,
        "sirolimus": True,
        "protac_1": True,
        "protac_2": True,
        "protac_3": False,
        "atorvastatin": False,
        "warfarin": False,
        "metoprolol": False
    }
    
    all_items = []
    with open("chameleon_local/benchmark.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] in targets:
                all_items.append((parts[0], parts[1]))
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] in targets:
                all_items.append((parts[0], parts[1]))

    print(f"{'Name':<15} | {'Label':<5} | {'Asph_ap':>7} | {'Asph_po':>7} | {'dAsph':>7} | {'Asph_min':>8} | {'Asph_max':>8}")
    print("-" * 80)

    for name, smi in all_items:
        try:
            mol = Chem.MolFromSmiles(smi)
            mol, cids = chameleon.embed_conformers(mol, 50, 4)
            energies = chameleon.minimize(mol, 4)
            keep, k_energies = chameleon.butina_prune(mol, energies, rms_cut=0.75)
            
            radii = chameleon.rdFreeSASA.classifyAtoms(mol)
            topo = chameleon.bond_path_distance_matrix(mol)
            
            psa_list, rg_list, imhb_list = [], [], []
            for cid in keep:
                psa_list.append(chameleon.compute_3d_psa(mol, cid, radii))
                rg_list.append(chameleon.compute_rg(mol, cid))
                imhb_list.append(chameleon.compute_imhb(mol, cid, topo))
            
            psa = np.asarray(psa_list)
            imhb = np.asarray(imhb_list, dtype=float)
            e_mmff = np.asarray(k_energies)
            
            e_apolar = e_mmff - chameleon.HBOND_EN_APOLAR * imhb
            e_polar = e_mmff - chameleon.POLAR_SOLV_KCAL_PER_A2 * psa + chameleon.HBOND_EN_POLAR * imhb
            w_apolar = chameleon.boltzmann_weights(e_apolar)
            w_polar = chameleon.boltzmann_weights(e_polar)
            
            asph_ap, asph_po, a_min, a_max, a_mean = get_shape_metrics(mol, keep, w_apolar, w_polar)
            dasph = asph_po - asph_ap
            
            is_pos = targets[name]
            l_str = "CHAM" if is_pos else "NON"
            print(f"{name:<15} | {l_str:<5} | {asph_ap:>7.3f} | {asph_po:>7.3f} | {dasph:>7.3f} | {a_min:>8.3f} | {a_max:>8.3f}")
        except Exception as e:
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    run_experiment()
