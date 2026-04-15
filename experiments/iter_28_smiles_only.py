import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

from rdkit import Chem
from rdkit.Chem import Descriptors
import json


def analyze_protac(smiles, name):
    """Analyze PROTAC from SMILES only (no 3D conformers)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    # Count atoms and bonds
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()

    # Count heteroatoms (N, O)
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    n_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_hetero = n_nitrogen + n_oxygen

    # Count specific bond types
    n_c_o_bonds = 0
    n_c_n_bonds = 0
    n_c_c_bonds = 0
    n_amide_bonds = 0
    n_ether_bonds = 0

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        n1 = a1.GetAtomicNum()
        n2 = a2.GetAtomicNum()

        # C-O bonds
        if (n1 == 6 and n2 == 8) or (n1 == 8 and n2 == 6):
            n_c_o_bonds += 1
            # Check for ether (C-O-C, not carbonyl)
            if not (a1.GetIsAromatic() or a2.GetIsAromatic()):
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    n_ether_bonds += 1

        # C-N bonds
        if (n1 == 6 and n2 == 7) or (n1 == 7 and n2 == 6):
            n_c_n_bonds += 1
            # Check for amide
            if bond.GetBondType() == Chem.BondType.SINGLE:
                # Check if adjacent to C=O
                for atom in [a1, a2]:
                    if atom.GetAtomicNum() == 6:  # Carbon
                        for neighbor in atom.GetNeighbors():
                            if (
                                neighbor.GetIdx() != a1.GetIdx()
                                and neighbor.GetIdx() != a2.GetIdx()
                            ):
                                bond2 = mol.GetBondBetweenAtoms(
                                    atom.GetIdx(), neighbor.GetIdx()
                                )
                                if (
                                    bond2
                                    and bond2.GetBondType() == Chem.BondType.DOUBLE
                                ):
                                    if neighbor.GetAtomicNum() == 8:  # C=O
                                        n_amide_bonds += 1

        # C-C bonds
        if n1 == 6 and n2 == 6:
            n_c_c_bonds += 1

    # Rotatable bonds
    n_rotatable = Descriptors.NumRotatableBonds(mol)

    # Ring information
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()

    # Molecular weight
    mw = Descriptors.MolWt(mol)

    return {
        "name": name,
        "mw": mw,
        "n_atoms": n_atoms,
        "n_bonds": n_bonds,
        "n_nitrogen": n_nitrogen,
        "n_oxygen": n_oxygen,
        "hetero_ratio": n_hetero / n_atoms if n_atoms > 0 else 0,
        "n_c_o": n_c_o_bonds,
        "n_c_n": n_c_n_bonds,
        "n_c_c": n_c_c_bonds,
        "n_ether": n_ether_bonds,
        "n_amide": n_amide_bonds,
        "ether_ratio": n_ether_bonds / n_bonds if n_bonds > 0 else 0,
        "rotatable": n_rotatable,
        "n_rings": n_rings,
        "rot_per_heavy": n_rotatable / n_atoms if n_atoms > 0 else 0,
        "flexibility_index": n_rotatable / n_rings if n_rings > 0 else n_rotatable,
    }


def main():
    print("=" * 70)
    print("ITERATION 28: SMILES-Based Topology Analysis (No 3D Conformers)")
    print("=" * 70)
    print()

    # Load PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) >= 2:
                name = parts[0]
                smiles = parts[1].strip()
                protacs.append((name, smiles))

    labels = {"protac_1": "CHAM", "protac_2": "CHAM", "protac_3": "NON"}

    results = []
    print("Processing PROTACs...")
    print("-" * 70)

    for name, smiles in protacs:
        result = analyze_protac(smiles, name)
        if result:
            results.append(result)
            label = labels.get(name, "?")
            print(f"\n{name} ({label}):")
            print(f"  MW: {result['mw']:.1f}")
            print(f"  Atoms: {result['n_atoms']}, Bonds: {result['n_bonds']}")
            print(f"  N: {result['n_nitrogen']}, O: {result['n_oxygen']}")
            print(f"  Hetero ratio: {result['hetero_ratio']:.3f}")
            print(
                f"  C-O bonds: {result['n_c_o']}, C-N: {result['n_c_n']}, C-C: {result['n_c_c']}"
            )
            print(
                f"  Ether bonds: {result['n_ether']} (ratio: {result['ether_ratio']:.3f})"
            )
            print(f"  Amide bonds: {result['n_amide']}")
            print(f"  Rotatable bonds: {result['rotatable']}")
            print(f"  Rings: {result['n_rings']}")
            print(f"  Flexibility index: {result['flexibility_index']:.2f}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    if len(results) >= 3:
        cham_ether = [
            r["ether_ratio"] for r in results if labels.get(r["name"]) == "CHAM"
        ]
        non_ether = [
            r["ether_ratio"] for r in results if labels.get(r["name"]) == "NON"
        ]

        if cham_ether and non_ether:
            cham_mean = sum(cham_ether) / len(cham_ether)
            non_mean = sum(non_ether) / len(non_ether)
            gap = cham_mean - non_mean

            print(f"\nChameleonic avg ether ratio: {cham_mean:.4f}")
            print(f"Non-chameleonic avg: {non_mean:.4f}")
            print(f"Separation: {gap:+.4f}")

            print("\nResults table:")
            print(
                f"{'Name':<12} {'Label':<6} {'Ether%':<8} {'RotBond':<8} {'FlexIdx':<8}"
            )
            print("-" * 50)
            for r in results:
                label = labels.get(r["name"], "?")
                print(
                    f"{r['name']:<12} {label:<6} {r['ether_ratio']:.3f}    {r['rotatable']:<8} {r['flexibility_index']:.2f}"
                )

            if gap > 0.01:
                print("\nVerdict: SUCCESS - PEG linkers have higher ether bond ratio")
            elif gap > 0:
                print("\nVerdict: PARTIAL SUCCESS - Correct ordering but small gap")
            else:
                print("\nVerdict: FAILED - Wrong ordering")

    # Write output
    with open("experiments/iter_28_output.txt", "w") as f:
        f.write("See stdout for results\n")

    # Write descriptors
    with open("experiments/iter_28_descriptors.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\nDescriptors written to experiments/iter_28_descriptors.json")


if __name__ == "__main__":
    main()
