import sys

# Force output
import os

os.environ["PYTHONUNBUFFERED"] = "1"

sys.path.insert(
    0, "C:/Users/mic23/prototype-embed/chameleon-research-loop/chameleon_local"
)

"""
Iteration 14: Unsatisfied H-Bond Donor Count

Hypothesis: Arvinas found unsatisfied HBD <= 2 is the most predictive rule for oral
PROTAC absorption (DOI: 10.1021/acs.jmedchem.3c00740). True chameleons mask polar surface
via IMHBs, leaving fewer HBDs exposed. We compute "unsatisfied HBD" = total HBD - IMHB_count
in the most compact conformer (minimum Rg).

Physical mechanism: In compact states, chameleons satisfy HBDs via IMHBs; non-chameleons
leave them exposed. This should be an inverse signal (lower unsatisfied HBD = more chameleonic).
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from chameleon import (
    embed_conformers,
    minimize,
    butina_prune,
    compute_3d_psa,
    compute_rg,
    compute_imhb,
    polar_atom_mask,
    bond_path_distance_matrix,
)


def compute_psa_rg_imhb_all_confs(mol, conf_ids):
    """Compute PSA, Rg, and IMHB for all conformers."""
    import numpy as np
    from rdkit.Chem import rdFreeSASA

    # Convert boolean mask to float for rdFreeSASA
    mask = polar_atom_mask(mol)
    radii = np.zeros(mol.GetNumAtoms(), dtype=np.float64)
    radii[mask] = 1.0  # dummy radius for polar atoms

    topo_dist = bond_path_distance_matrix(mol)
    results = []
    for cid in conf_ids:
        # CalcSASA and sum polar contributions
        rdFreeSASA.CalcSASA(mol, radii, confIdx=cid)
        pmask = polar_atom_mask(mol)
        psa = 0.0
        for atom in mol.GetAtoms():
            if pmask[atom.GetIdx()]:
                psa += float(atom.GetProp("SASA"))

        rg = compute_rg(mol, cid)
        imhb = compute_imhb(mol, cid, topo_dist)
        results.append({"psa": psa, "rg": rg, "imhb": imhb})
    return results


def count_hbd(smiles: str) -> int:
    """Count H-bond donors using Lipinski definition."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    mol = Chem.AddHs(mol)
    return rdMolDescriptors.CalcNumHBD(mol)


def analyze_molecule(smiles: str, name: str, n_conf: int = 50):
    """Compute unsatisfied HBD for a molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Generate conformers
    mol, cids = embed_conformers(mol, n_conf=n_conf, n_threads=4)
    energies = minimize(mol, n_threads=4)

    # Prune
    keep_ids, keep_energies = butina_prune(mol, energies, rms_cut=0.75)

    if len(keep_ids) == 0:
        return None

    # Compute descriptors
    results = compute_psa_rg_imhb_all_confs(mol, keep_ids)

    # Find most compact conformer (minimum Rg)
    rg_values = [r["rg"] for r in results]
    min_rg_idx = np.argmin(rg_values)

    # Get IMHB count in most compact conformer
    imhb_compact = results[min_rg_idx]["imhb"]

    # Total HBD
    total_hbd = count_hbd(smiles)

    # Unsatisfied HBD
    unsatisfied_hbd = max(0, total_hbd - imhb_compact)

    return {
        "name": name,
        "smiles": smiles,
        "total_hbd": total_hbd,
        "imhb_compact": imhb_compact,
        "unsatisfied_hbd": unsatisfied_hbd,
        "min_rg": results[min_rg_idx]["rg"],
        "max_rg": max(rg_values),
        "rg_ratio": max(rg_values) / min(rg_values),
    }


def main():
    output = []

    def log(msg):
        output.append(msg)
        print(msg)

    # Test PROTACs
    protacs = [
        ("COc1cc2c(cc1OC)C(=O)N(CCCN1CCC(C(=O)Nc3ccccc3)CC1)C2", "protac_1", "CHAM"),
        ("COc1cc2c(cc1OC)C(=O)N(CCCCN1CCC(C(=O)Nc3ccccc3)CC1)C2", "protac_2", "CHAM"),
        ("COc1cc2c(cc1OC)C(=O)N(CCCCCC(=O)Nc1ccccc1)C2", "protac_3", "NON"),
    ]

    log("=" * 70)
    log("Iteration 14: Unsatisfied H-Bond Donor Count")
    log("=" * 70)
    log(f"\nN_CONF = 50 per molecule")
    log(f"\nHypothesis: Lower unsatisfied HBD in compact state = more chameleonic")
    log(f"(Arvinas threshold: unsatisfied HBD <= 2 for oral absorption)")
    log("")

    results = []
    for smiles, name, label in protacs:
        result = analyze_molecule(smiles, name, n_conf=50)
        if result:
            result["label"] = label
            results.append(result)

    # Print results
    log(
        f"\n{'Name':<15} {'Label':<6} {'Total HBD':<10} {'IMHB_compact':<12} {'Unsat HBD':<10} {'Rg_min':<8} {'Rg_max':<8}"
    )
    log("-" * 85)

    for r in results:
        log(
            f"{r['name']:<15} {r['label']:<6} {r['total_hbd']:<10} {r['imhb_compact']:<12.2f} {r['unsatisfied_hbd']:<10.2f} {r['min_rg']:<8.2f} {r['max_rg']:<8.2f}"
        )

    # Calculate signal strength
    log("\n" + "=" * 70)
    log("Signal Analysis:")
    log("=" * 70)

    cham_unsat = [r["unsatisfied_hbd"] for r in results if r["label"] == "CHAM"]
    non_unsat = [r["unsatisfied_hbd"] for r in results if r["label"] == "NON"]

    if cham_unsat and non_unsat:
        cham_mean = np.mean(cham_unsat)
        non_mean = np.mean(non_unsat)

        log(f"\nMean unsatisfied HBD (Chameleonic): {cham_mean:.2f}")
        log(f"Mean unsatisfied HBD (Non-chameleonic): {non_mean:.2f}")
        log(f"Difference: {non_mean - cham_mean:.2f}")

        # Arvinas threshold test
        log(f"\n--- Arvinas Threshold Test (unsatisfied HBD <= 2) ---")
        for r in results:
            status = (
                "PASS"
                if (r["label"] == "CHAM" and r["unsatisfied_hbd"] <= 2)
                or (r["label"] == "NON" and r["unsatisfied_hbd"] > 2)
                else "FAIL"
            )
            log(
                f"{r['name']:<15}: unsat_HBD={r['unsatisfied_hbd']:.2f}, label={r['label']}, threshold_pass={r['unsatisfied_hbd'] <= 2} [{status}]"
            )

        # Check if unsatisfied HBD discriminates protac_3
        protac_3_unsat = [
            r["unsatisfied_hbd"] for r in results if r["name"] == "protac_3"
        ][0]
        protac_1_unsat = [
            r["unsatisfied_hbd"] for r in results if r["name"] == "protac_1"
        ][0]
        protac_2_unsat = [
            r["unsatisfied_hbd"] for r in results if r["name"] == "protac_2"
        ][0]

        log(f"\n--- PROTAC Separation ---")
        log(f"protac_1 (CHAM): {protac_1_unsat:.2f}")
        log(f"protac_2 (CHAM): {protac_2_unsat:.2f}")
        log(f"protac_3 (NON):  {protac_3_unsat:.2f}")

        if protac_3_unsat > max(protac_1_unsat, protac_2_unsat):
            log(
                "\nSUCCESS: protac_3 has HIGHER unsatisfied HBD than chameleonic PROTACs"
            )
            log("  (Non-chameleons leave more HBDs exposed)")
        elif protac_3_unsat < min(protac_1_unsat, protac_2_unsat):
            log("\nFAIL: protac_3 has LOWER unsatisfied HBD than chameleonic PROTACs")
        else:
            log("\nUNCLEAR: protac_3 unsatisfied HBD overlaps with chameleonic range")

    log("\n" + "=" * 70)
    log("Summary:")
    log("=" * 70)
    log("Unsatisfied HBD = Total HBD - IMHB in most compact conformer")
    log("Lower values indicate more HBDs are satisfied via IMHBs (chameleonic)")
    log("")
    log("Arvinas finding for oral absorption: unsatisfied HBD <= 2")
    log("Expected for chameleonicity: similar pattern - chameleons mask HBDs")

    # Write results to file for capture
    with open("experiments/iter_14_output.txt", "w") as f:
        f.write("\n".join(output))


if __name__ == "__main__":
    main()
