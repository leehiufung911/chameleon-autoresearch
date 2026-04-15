
import sys
import os
import math
import numpy as np
import time
import traceback
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdFreeSASA, Descriptors

# --- Manual Butina implementation to avoid failing rdkit.ML.Cluster import ---
def manual_butina_cluster_data(data, npts, threshold, isDistData=True, reordering=True):
    neighbors = [[] for _ in range(npts)]
    idx = 0
    for i in range(1, npts):
        for j in range(i):
            d = data[idx]
            if d <= threshold:
                neighbors[i].append(j)
                neighbors[j].append(i)
            idx += 1
    order = sorted(range(npts), key=lambda i: len(neighbors[i]), reverse=True)
    clusters = []
    assigned = set()
    for i in order:
        if i not in assigned:
            cluster = [i]
            assigned.add(i)
            for j in neighbors[i]:
                if j not in assigned:
                    cluster.append(j)
                    assigned.add(j)
            clusters.append(tuple(cluster))
    return tuple(clusters)

# Mock the failing module
from unittest.mock import MagicMock
mock_butina = MagicMock()
mock_butina.ClusterData = manual_butina_cluster_data
sys.modules["rdkit.ML.Cluster"] = MagicMock()
sys.modules["rdkit.ML.Cluster.Butina"] = mock_butina

# Now we can import from chameleon
sys.path.append(os.getcwd())
from chameleon_local.chameleon import (
    embed_conformers, minimize, butina_prune, 
    compute_3d_psa, compute_rg, compute_imhb, 
    boltzmann_weights, bond_path_distance_matrix,
    HBOND_EN_APOLAR, HBOND_EN_POLAR, POLAR_SOLV_KCAL_PER_A2
)

def compute_shape_descriptors(mol, conf_id):
    """Compute various shape descriptors for a given conformer."""
    asphericity = rdMolDescriptors.CalcAsphericity(mol, confId=conf_id)
    eccentricity = rdMolDescriptors.CalcEccentricity(mol, confId=conf_id)
    sphericality = rdMolDescriptors.CalcSphericalityIndex(mol, confId=conf_id)
    return asphericity, eccentricity, sphericality

def run_experiment(name, smiles, n_conf=100):
    try:
        mol0 = Chem.MolFromSmiles(smiles)
        if mol0 is None: return None
        mw = Descriptors.MolWt(mol0)
        
        mol, cids = embed_conformers(mol0, n_conf, n_threads=0)
        energies = minimize(mol, n_threads=0)
        keep, k_energies = butina_prune(mol, energies, rms_cut=0.75)
        
        radii = rdFreeSASA.classifyAtoms(mol)
        topo = bond_path_distance_matrix(mol)
        
        psa_list, rg_list, imhb_list = [], [], []
        asph_list, ecc_list, sph_list = [], [], []
        
        for cid in keep:
            psa_list.append(compute_3d_psa(mol, cid, radii))
            rg_list.append(compute_rg(mol, cid))
            imhb_list.append(compute_imhb(mol, cid, topo))
            asph, ecc, sph = compute_shape_descriptors(mol, cid)
            asph_list.append(asph)
            ecc_list.append(ecc)
            sph_list.append(sph)
            
        psa = np.array(psa_list)
        rg = np.array(rg_list)
        imhb = np.array(imhb_list, dtype=float)
        e_mmff = np.array(k_energies)
        asph = np.array(asph_list)
        sph = np.array(sph_list)
        
        # Ensembles
        e_apolar = e_mmff - HBOND_EN_APOLAR * imhb
        e_polar = e_mmff - POLAR_SOLV_KCAL_PER_A2 * psa + HBOND_EN_POLAR * imhb
        w_apolar = boltzmann_weights(e_apolar)
        w_polar = boltzmann_weights(e_polar)
        
        # Delta Shape
        sph_ap = np.sum(w_apolar * sph)
        sph_po = np.sum(w_polar * sph)
        dSph = sph_ap - sph_po  # Positive means more spherical in apolar (folded)
        
        asph_ap = np.sum(w_apolar * asph)
        asph_po = np.sum(w_polar * asph)
        dAsph = asph_po - asph_ap # Positive means more asymmetric in polar (extended)
        
        # CI original
        psa_range = psa.max() - psa.min()
        rg_ratio = rg.max() / max(rg.min(), 1e-6)
        term_psa = psa_range / math.sqrt(mw)
        term_imhb = imhb.mean()
        term_rg = math.log(rg_ratio)
        ci_orig = 2.0 * term_psa + 1.0 * term_imhb + 3.0 * term_rg
        
        return {
            "name": name,
            "ci_orig": ci_orig,
            "sph_ap": sph_ap,
            "sph_po": sph_po,
            "dSph": dSph,
            "asph_ap": asph_ap,
            "asph_po": asph_po,
            "dAsph": dAsph,
            "sph_mean": sph.mean(),
            "asph_mean": asph.mean()
        }
    except Exception as e:
        print(f"Error in {name}: {e}")
        traceback.print_exc()
        return None

if __name__ == "__main__":
    test_mols = [
        ("protac_1", "O=C(C(N1C(C2=CC=CC(NCCOCCOCCOCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=C2C1=O)=O)CC7)NC7=O", "Cham"),
        ("protac_2", "O=C(C(N1C(C2=CC=CC(NC(COCCOCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O", "Cham"),
        ("protac_3", "O=C(C(N1C(C2=CC=CC(OCC(NCCCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O", "Non"),
        ("cyclosporinA", "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C", "Cham"),
        ("atorvastatin", "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4", "Non")
    ]
    
    print(f"{'Name':<15} | {'GT':<5} | {'CI':>5} | {'dSph':>7} | {'dAsph':>7} | {'Sph_m':>7} | {'Asph_m':>7}")
    print("-" * 75)
    
    results = []
    for name, smi, gt in test_mols:
        res = run_experiment(name, smi, n_conf=100)
        if res:
            res["gt"] = gt
            results.append(res)
            print(f"{name:<15} | {gt:<5} | {res['ci_orig']:5.2f} | {res['dSph']:7.4f} | {res['dAsph']:7.4f} | {res['sph_mean']:7.4f} | {res['asph_mean']:7.4f}")
