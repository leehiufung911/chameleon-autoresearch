"""
Iteration 15: Combined Veto Pipeline

Hypothesis: Merge the best findings from previous iterations into a single
pipeline that applies structural vetoes based on physical reasoning, combined
with multiplicative CI coupling and saturation normalization.

Components:
1. MaxCPath >= 5 + IMHB_mean < 1.0 → alkyl veto (iteration 13)
2. Multiplicative CI coupling: (dPSA + dRg) * IMHB_saturation (iteration 6)
3. IMHB saturation normalization: IMHB_mean / N_hbond_sites (iteration 9)

Physical Mechanism: True chameleonicity requires both flexible folding machinery
(multiplicative coupling of dPSA/dRg with IMHB saturation) AND a structural
bias against long alkyl chains that prefer extended states in apolar solvent.
"""

import sys
import os
import math

# Add chameleon_local to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "chameleon_local"))

# Test imports
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem, rdFreeSASA, rdMolDescriptors
    from rdkit.Chem import rdMolTransforms
    from rdkit.ML.Cluster import Butina
    import numpy as np

    RDKIT_AVAILABLE = True
except ImportError as e:
    print(f"Import error: {e}")
    RDKIT_AVAILABLE = False
    Chem = None
    Descriptors = None
    rdMolDescriptors = None


def compute_max_contiguous_carbon_path(mol):
    """Find longest contiguous non-ring carbon chain."""
    if mol is None or Chem is None:
        return 0

    ring_atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms.update(ring)

    carbon_atoms = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "C" and atom.GetIdx() not in ring_atoms
    ]

    if not carbon_atoms:
        return 0

    # Build adjacency for non-ring carbons
    adj = {idx: [] for idx in carbon_atoms}
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a1 in adj and a2 in adj:
            adj[a1].append(a2)
            adj[a2].append(a1)

    # DFS to find longest path
    def dfs(node, visited):
        visited.add(node)
        max_len = 1
        for neighbor in adj[node]:
            if neighbor not in visited:
                max_len = max(max_len, 1 + dfs(neighbor, visited))
        visited.remove(node)
        return max_len

    max_path = 0
    for start in carbon_atoms:
        max_path = max(max_path, dfs(start, set()))

    return max_path


def count_hbond_sites(mol):
    """Count potential H-bond donors and acceptors (topo_dist >= 4)."""
    if mol is None or Chem is None:
        return 1

    donors = []
    acceptors = []

    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in ("N", "O"):
            acceptors.append(atom.GetIdx())
            for n in atom.GetNeighbors():
                if n.GetSymbol() == "H":
                    donors.append((atom.GetIdx(), n.GetIdx()))

    topo = Chem.GetDistanceMatrix(mol)
    count = 0
    for d_idx, h_idx in donors:
        for a_idx in acceptors:
            if a_idx == d_idx:
                continue
            if topo[d_idx, a_idx] >= 4:
                count += 1

    return max(count, 1)  # Avoid division by zero


def compute_imhb_for_mol(mol):
    """Compute IMHB for a molecule (simplified)."""
    if mol is None or Chem is None:
        return 0

    topo = Chem.GetDistanceMatrix(mol)
    donor_pairs = []
    acceptors = []

    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in ("N", "O"):
            acceptors.append(atom.GetIdx())
            for n in atom.GetNeighbors():
                if n.GetSymbol() == "H":
                    donor_pairs.append((atom.GetIdx(), n.GetIdx()))

    if not donor_pairs or not acceptors:
        return 0

    # Simplified: count potential H-bond sites as proxy
    return min(len(donor_pairs), len(acceptors)) // 2


def estimate_ci_components(smiles, name):
    """Estimate CI components from 2D structure."""
    if Chem is None or Descriptors is None or rdMolDescriptors is None:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    mw = Descriptors.MolWt(mol)

    # Structural descriptors
    max_c_path = compute_max_contiguous_carbon_path(mol)
    n_sites = count_hbond_sites(mol)
    imhb_estimate = compute_imhb_for_mol(mol)

    # Estimate IMHB mean from 2D (rough proxy)
    imhb_mean = imhb_estimate / max(1, n_sites / 10)  # Normalized estimate

    # 2D TPSA as proxy for PSA range
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    psa_range = tpsa * 0.3  # Estimated range from 2D

    # Heavy atoms for Rg estimate
    heavy = mol.GetNumHeavyAtoms()
    rg_ratio = 1 + (heavy / 50)  # Estimated expansion ratio

    # CI components
    term_psa = psa_range / math.sqrt(max(mw, 1.0))
    term_rg = math.log(max(rg_ratio, 1.0))

    # IMHB saturation
    imhb_saturation = imhb_mean / n_sites if n_sites > 0 else 0

    # Original CI
    ci_orig = 2.0 * term_psa + 1.0 * imhb_mean + 3.0 * term_rg

    # Coupled CI with saturation
    ci_coupled = (term_psa + 2.0 * term_rg) * (1.0 + 10.0 * imhb_saturation)

    # Apply alkyl linker veto (MaxCPath >= 5 AND IMHB_mean < 1.0)
    veto_applied = max_c_path >= 5 and imhb_mean < 1.0
    ci_final = ci_coupled * 0.5 if veto_applied else ci_coupled

    return {
        "name": name,
        "smiles": smiles,
        "mw": mw,
        "max_c_path": max_c_path,
        "n_sites": n_sites,
        "imhb_mean": imhb_mean,
        "imhb_saturation": imhb_saturation,
        "term_psa": term_psa,
        "term_rg": term_rg,
        "ci_original": ci_orig,
        "ci_coupled": ci_coupled,
        "ci_final": ci_final,
        "veto_applied": veto_applied,
    }


def main():
    print("=" * 80)
    print("Iteration 15: Combined Veto Pipeline (2D-based estimation)")
    print("=" * 80)
    print("\nNOTE: Using 2D descriptors as fast proxies for 3D conformer analysis")
    print("Full 3D pipeline requires MMFF minimization which is slow.")
    print("This tests the veto logic and formula structure.")
    print()

    # Test on the three PROTACs
    protacs = [
        (
            "protac_1",
            "O=C(C(N1C(C2=CC=CC(NCCOCCOCCOCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=C2C1=O)=O)CC7)NC7=O",
        ),
        (
            "protac_2",
            "O=C(C(N1C(C2=CC=CC(NC(COCCOCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
        ),
        (
            "protac_3",
            "O=C(C(N1C(C2=CC=CC(OCC(NCCCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
        ),
    ]

    results = []
    for name, smiles in protacs:
        result = estimate_ci_components(smiles, name)
        if result:
            results.append(result)

            print(f"\n{name}:")
            print(f"  MW: {result['mw']:.1f}")
            print(
                f"  MaxCPath: {result['max_c_path']}, IMHB_mean: {result['imhb_mean']:.2f}"
            )
            print(
                f"  IMHB_saturation: {result['imhb_saturation']:.4f} (n_sites: {result['n_sites']})"
            )
            print(f"  CI_original: {result['ci_original']:.2f}")
            print(f"  CI_coupled:  {result['ci_coupled']:.2f}")
            print(f"  Veto_applied: {result['veto_applied']}")
            print(f"  CI_final:    {result['ci_final']:.2f}")

            # Verdict
            if result["ci_final"] >= 4.0:
                verdict = "CHAMELEONIC"
            elif result["ci_final"] >= 2.0:
                verdict = "MARGINAL"
            else:
                verdict = "NON-CHAMELEONIC"
            print(f"  -> {verdict}")

    print("\n" + "=" * 80)
    print("Summary:")
    print("=" * 80)
    print(
        f"{'Name':<12} {'Original':>10} {'Coupled':>10} {'Final':>10} {'Veto':>6} {'Verdict':>15}"
    )
    print("-" * 80)

    ground_truth = {"protac_1": "CHAM", "protac_2": "CHAM", "protac_3": "NON"}

    for r in results:
        gt = ground_truth.get(r["name"], "UNK")
        if r["ci_final"] >= 4.0:
            pred = "CHAM"
        elif r["ci_final"] >= 2.0:
            pred = "MARG"
        else:
            pred = "NON"

        correct = "✓" if pred == gt else "✗"
        print(
            f"{r['name']:<12} {r['ci_original']:>10.2f} {r['ci_coupled']:>10.2f} "
            f"{r['ci_final']:>10.2f} {str(r['veto_applied']):>6} {pred:>14} {correct}"
        )

    print("\nExpected: protac_1=CHAM, protac_2=CHAM, protac_3=NON")
    print("\nKey Finding:")

    # Check if protac_3 was correctly vetoed
    p3 = next((r for r in results if r["name"] == "protac_3"), None)
    if p3 and p3["veto_applied"]:
        print(
            f"  ✓ Veto correctly applied to protac_3 (MaxCPath={p3['max_c_path']}, IMHB={p3['imhb_mean']:.2f})"
        )
        if p3["ci_final"] < 4.0:
            print(
                f"  ✓ protac_3 classified as NON-CHAMELEONIC (CI={p3['ci_final']:.2f})"
            )
        else:
            print(f"  ⚠ Veto applied but CI still high ({p3['ci_final']:.2f})")

    # Check if chameleons preserved
    p1 = next((r for r in results if r["name"] == "protac_1"), None)
    p2 = next((r for r in results if r["name"] == "protac_2"), None)
    if p1 and not p1["veto_applied"] and p1["ci_final"] >= 4.0:
        print(f"  ✓ protac_1 remains CHAMELEONIC (CI={p1['ci_final']:.2f})")
    if p2 and not p2["veto_applied"] and p2["ci_final"] >= 4.0:
        print(f"  ✓ protac_2 remains CHAMELEONIC (CI={p2['ci_final']:.2f})")


if __name__ == "__main__":
    if not RDKIT_AVAILABLE:
        print("RDKit not available - this experiment requires RDKit")
        sys.exit(1)
    main()
