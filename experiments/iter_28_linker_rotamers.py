import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
import json
from collections import defaultdict


def get_linker_rotatable_bonds(mol, linker_atoms):
    """
    Identify rotatable bonds within the linker region.
    Returns list of (atom1_idx, atom2_idx, bond_type) tuples.
    """
    linker_set = set(linker_atoms)
    rotatable_bonds = []

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()

        # Both atoms must be in linker
        if a1 not in linker_set or a2 not in linker_set:
            continue

        # Skip non-rotatable bonds (ring, double, triple, amide)
        if bond.IsInRing():
            continue
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        # Check for amide (C=O next to N)
        atom1 = mol.GetAtomWithIdx(a1)
        atom2 = mol.GetAtomWithIdx(a2)

        if is_amide_bond(mol, a1, a2):
            continue

        # Classify bond type based on atoms
        bond_type = classify_bond_type(atom1, atom2)
        rotatable_bonds.append((a1, a2, bond_type))

    return rotatable_bonds


def is_amide_bond(mol, idx1, idx2):
    """Check if bond is an amide C-N bond."""
    atom1 = mol.GetAtomWithIdx(idx1)
    atom2 = mol.GetAtomWithIdx(idx2)

    # Check if C-N bond with adjacent C=O
    if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7) or (
        atom1.GetAtomicNum() == 7 and atom2.GetAtomicNum() == 6
    ):
        # Check for adjacent carbonyl
        c_idx = idx1 if atom1.GetAtomicNum() == 6 else idx2
        c_atom = mol.GetAtomWithIdx(c_idx)

        for neighbor in c_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx == idx1 or n_idx == idx2:
                continue
            bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                if neighbor.GetAtomicNum() == 8:  # C=O
                    return True
    return False


def classify_bond_type(atom1, atom2):
    """
    Classify bond based on atom types.
    Returns: 'CC', 'CO', 'OC', 'CN', 'NC', or 'CX' (generic)
    """
    n1 = atom1.GetAtomicNum()
    n2 = atom2.GetAtomicNum()

    # Carbon-carbon
    if n1 == 6 and n2 == 6:
        return "CC"

    # Carbon-oxygen (either direction)
    if (n1 == 6 and n2 == 8) or (n1 == 8 and n2 == 6):
        return "CO"

    # Carbon-nitrogen
    if (n1 == 6 and n2 == 7) or (n1 == 7 and n2 == 6):
        return "CN"

    return "CX"


def measure_dihedral_angle(mol, a1, a2):
    """
    Measure the dihedral angle around a rotatable bond.
    Need 4 atoms: neighbor(a1) - a1 - a2 - neighbor(a2)
    """
    atom1 = mol.GetAtomWithIdx(a1)
    atom2 = mol.GetAtomWithIdx(a2)

    # Find neighbors not part of the bond
    neighbors1 = [n.GetIdx() for n in atom1.GetNeighbors() if n.GetIdx() != a2]
    neighbors2 = [n.GetIdx() for n in atom2.GetNeighbors() if n.GetIdx() != a1]

    if not neighbors1 or not neighbors2:
        return None

    # Take first neighbor
    prev = neighbors1[0]
    next_a = neighbors2[0]

    try:
        conf = mol.GetConformer()
        angle = rdMolDescriptors.CalcDihedral(conf, prev, a1, a2, next_a)
        return angle  # in radians
    except:
        return None


def classify_rotamer(angle_rad, bond_type):
    """
    Classify dihedral angle into rotamer state.
    Uses standard rotamer definitions from polymer physics:
    - anti: 180 deg +/- 30 deg
    - gauche+: 60 deg +/- 30 deg
    - gauche-: -60 deg +/- 30 deg
    """
    if angle_rad is None:
        return "unknown"

    angle_deg = np.degrees(angle_rad)

    # Normalize to -180 to 180
    while angle_deg > 180:
        angle_deg -= 360
    while angle_deg < -180:
        angle_deg += 360

    if 150 <= abs(angle_deg) <= 180:
        return "anti"
    elif 30 <= angle_deg <= 90:
        return "gauche_plus"
    elif -90 <= angle_deg <= -30:
        return "gauche_minus"
    else:
        return "other"


def identify_warheads(mol):
    """
    Identify warhead/ligand regions (largest ring systems).
    Returns list of atom indices in warhead regions.
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    if not atom_rings:
        # No rings - use aromatic atoms
        return [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()]

    # Find largest rings (likely warheads)
    ring_sizes = [(len(r), r) for r in atom_rings]
    ring_sizes.sort(reverse=True)

    # Take top 2 largest rings as warheads
    warhead_atoms = set()
    for size, ring in ring_sizes[:2]:
        warhead_atoms.update(ring)

    return list(warhead_atoms)


def identify_linker_region(mol):
    """
    Identify linker atoms between warheads.
    Simplified: atoms not in the two largest rings.
    """
    warhead_atoms = set(identify_warheads(mol))

    # Linker is everything else
    linker_atoms = [i for i in range(mol.GetNumAtoms()) if i not in warhead_atoms]

    return linker_atoms


def compute_rotamer_distribution(smiles, name, n_conformers=30):
    """
    Compute rotamer distribution for linker bonds.
    This simulates what would be extracted from PDB structures.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    # Generate diverse conformers
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    ps.numThreads = 0

    AllChem.EmbedMultipleConfs(mol, n_conformers, ps)

    if mol.GetNumConformers() == 0:
        return None

    # Minimize each conformer
    for conf_id in range(mol.GetNumConformers()):
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, mmffVariant="MMFF94s")
        except:
            pass

    # Identify linker region
    linker_atoms = identify_linker_region(mol)

    if len(linker_atoms) < 3:
        return None

    # Get rotatable bonds in linker
    rotatable_bonds = get_linker_rotatable_bonds(mol, linker_atoms)

    if not rotatable_bonds:
        return None

    # Count rotamers across conformers
    rotamer_counts = defaultdict(lambda: defaultdict(int))

    for conf_id in range(mol.GetNumConformers()):
        mol.SetActiveConf(conf_id)

        for a1, a2, bond_type in rotatable_bonds:
            angle = measure_dihedral_angle(mol, a1, a2)
            rotamer = classify_rotamer(angle, bond_type)

            key = (a1, a2, bond_type)
            rotamer_counts[key][rotamer] += 1

    # Compute rotamer preferences
    total_bonds = len(rotatable_bonds)

    # Aggregate by bond type
    type_counts = defaultdict(lambda: defaultdict(int))
    type_totals = defaultdict(int)

    for key, counts in rotamer_counts.items():
        bond_type = key[2]
        for rotamer, count in counts.items():
            type_counts[bond_type][rotamer] += count
            type_totals[bond_type] += count

    # Calculate preferences per bond type
    preferences = {}
    for bond_type in ["CC", "CO", "CN"]:
        if type_totals[bond_type] > 0:
            anti_frac = type_counts[bond_type].get("anti", 0) / type_totals[bond_type]
            gauche_frac = (
                type_counts[bond_type].get("gauche_plus", 0)
                + type_counts[bond_type].get("gauche_minus", 0)
            ) / type_totals[bond_type]

            preferences[bond_type] = {
                "anti_fraction": anti_frac,
                "gauche_fraction": gauche_frac,
                "gauche_preference": gauche_frac
                - anti_frac,  # Positive = prefers gauche
                "total_observations": type_totals[bond_type],
            }

    # Compute linker-level metrics
    cc_pref = preferences.get("CC", {}).get("gauche_preference", 0)
    co_pref = preferences.get("CO", {}).get("gauche_preference", 0)

    # Weighted average by bond count
    cc_count = type_totals.get("CC", 0)
    co_count = type_totals.get("CO", 0)
    total_count = cc_count + co_count

    if total_count > 0:
        overall_gauche_preference = (
            cc_pref * cc_count + co_pref * co_count
        ) / total_count
    else:
        overall_gauche_preference = 0

    return {
        "name": name,
        "n_conformers": mol.GetNumConformers(),
        "n_rotatable_bonds": total_bonds,
        "linker_atoms": len(linker_atoms),
        "preferences_by_type": preferences,
        "overall_gauche_preference": overall_gauche_preference,
        "cc_bonds": cc_count,
        "co_bonds": co_count,
        "cn_bonds": type_totals.get("CN", 0),
        "ether_ratio": co_count / total_count if total_count > 0 else 0,
        "rotamer_entropy": compute_rotamer_entropy(rotamer_counts)
        if rotamer_counts
        else 0,
    }


def compute_rotamer_entropy(rotamer_counts):
    """Compute Shannon entropy of rotamer distribution."""
    total = sum(sum(counts.values()) for counts in rotamer_counts.values())
    if total == 0:
        return 0

    entropy = 0
    for key, counts in rotamer_counts.items():
        bond_total = sum(counts.values())
        for rotamer, count in counts.items():
            p = count / bond_total
            if p > 0:
                entropy -= p * np.log2(p) * (bond_total / total)  # Weight by bond

    return entropy


def main():
    output_lines = [
        "=" * 70,
        "ITERATION 28: Rotamer Library Analysis from Conformer Ensembles",
        "=" * 70,
        "",
        "Hypothesis: PEG linkers (gauche preference) vs alkyl linkers (anti preference)",
        "can be detected by analyzing dihedral angle distributions. This simulates",
        "extracting empirical rotamer distributions from PDB PROTAC structures.",
        "",
        "Physical Mechanism:",
        "- Ether bonds (C-O) favor gauche (~60 deg, -60 deg) due to anomeric effect",
        "- Carbon bonds (C-C) in alkanes favor anti (~180 deg) from sterics",
        "- Chameleonic PROTACs have more ether linkers -> higher gauche preference",
        "- Non-chameleonic alkyl linkers have more C-C bonds -> anti preference",
        "",
    ]

    # Load PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) >= 2:
                name = parts[0]
                smiles = parts[1].strip()
                protacs.append((name, smiles))

    # Labels (from paper)
    labels = {"protac_1": "CHAM", "protac_2": "CHAM", "protac_3": "NON"}

    results = []

    output_lines.append("Processing PROTACs...")
    output_lines.append("-" * 70)

    for name, smiles in protacs:
        result = compute_rotamer_distribution(smiles, name, n_conformers=30)

        if result:
            results.append(result)

            output_lines.append(f"\n{name} ({labels.get(name, 'unknown')}):")
            output_lines.append(f"  Linker atoms: {result['linker_atoms']}")
            output_lines.append(f"  Rotatable bonds: {result['n_rotatable_bonds']}")
            output_lines.append(
                f"  CC bonds: {result['cc_bonds']}, CO bonds: {result['co_bonds']}"
            )
            output_lines.append(f"  Ether ratio: {result['ether_ratio']:.3f}")
            output_lines.append(
                f"  Overall gauche preference: {result['overall_gauche_preference']:.3f}"
            )
            output_lines.append(f"  Rotamer entropy: {result['rotamer_entropy']:.3f}")

            # Bond type details
            for bond_type, prefs in result["preferences_by_type"].items():
                output_lines.append(
                    f"  {bond_type}: anti={prefs['anti_fraction']:.2f}, "
                    f"gauche={prefs['gauche_fraction']:.2f}, "
                    f"pref={prefs['gauche_preference']:+.3f}"
                )
        else:
            output_lines.append(f"\n{name}: FAILED to analyze")

    # Summary
    output_lines.append("\n" + "=" * 70)
    output_lines.append("SUMMARY")
    output_lines.append("=" * 70)

    # Extract key metrics
    cham_values = []
    non_values = []

    for r in results:
        if labels.get(r["name"]) == "CHAM":
            cham_values.append(r["overall_gauche_preference"])
        else:
            non_values.append(r["overall_gauche_preference"])

    if cham_values and non_values:
        cham_mean = np.mean(cham_values)
        non_mean = np.mean(non_values)
        gap = cham_mean - non_mean

        output_lines.append(f"\nChameleonic avg gauche preference: {cham_mean:.3f}")
        output_lines.append(f"Non-chameleonic avg: {non_mean:.3f}")
        output_lines.append(f"Separation gap: {gap:+.3f} (positive = correct)")

        # Check individual results
        output_lines.append("\nIndividual results:")
        for r in results:
            label = labels.get(r["name"], "?")
            pref = r["overall_gauche_preference"]
            output_lines.append(
                f"  {r['name']:<12} ({label}): gauche_pref = {pref:+.3f}"
            )

        if gap > 0:
            output_lines.append(
                "\nVerdict: SUCCESS - PEG linkers show higher gauche preference"
            )
        else:
            output_lines.append(
                "\nVerdict: FAILED - No separation between linker types"
            )

    # Write descriptors JSON for benchmark
    descriptors = []
    for r in results:
        descriptors.append(
            {
                "name": r["name"],
                "label": labels.get(r["name"], "unknown"),
                "overall_gauche_preference": r["overall_gauche_preference"],
                "ether_ratio": r["ether_ratio"],
                "rotamer_entropy": r["rotamer_entropy"],
                "cc_bonds": r["cc_bonds"],
                "co_bonds": r["co_bonds"],
                "n_rotatable_bonds": r["n_rotatable_bonds"],
            }
        )

    with open("experiments/iter_28_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)

    final_output = "\n".join(output_lines)

    # Write output file
    with open("experiments/iter_28_output.txt", "w") as f:
        f.write(final_output)

    print(final_output)
    print("\nOutput written to experiments/iter_28_output.txt")
    print("Descriptors written to experiments/iter_28_descriptors.json")


if __name__ == "__main__":
    main()
