import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

print("=" * 70)
print("ITERATION 28: Minimal SMILES Analysis")
print("=" * 70)
print()

from rdkit import Chem
from rdkit.Chem import Descriptors
import json

# Load first PROTAC
with open("chameleon_local/user_protacs.tsv", "r") as f:
    line = f.readline()
    parts = line.strip().split(maxsplit=1)
    name = parts[0]
    smiles = parts[1].strip()

print(f"Processing {name}")
print(f"SMILES: {smiles[:50]}...")

mol = Chem.MolFromSmiles(smiles)
if mol:
    mol = Chem.AddHs(mol)
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()
    mw = Descriptors.MolWt(mol)
    rot = Descriptors.NumRotatableBonds(mol)

    # Count C-O bonds
    n_co = 0
    for bond in mol.GetBonds():
        n1 = bond.GetBeginAtom().GetAtomicNum()
        n2 = bond.GetEndAtom().GetAtomicNum()
        if (n1 == 6 and n2 == 8) or (n1 == 8 and n2 == 6):
            n_co += 1

    print(f"\nResults:")
    print(f"  Atoms: {n_atoms}, Bonds: {n_bonds}")
    print(f"  MW: {mw:.1f}")
    print(f"  Rotatable bonds: {rot}")
    print(f"  C-O bonds: {n_co}")
    print(f"  C-O ratio: {n_co / n_bonds:.3f}")

    # Save results
    result = {
        "name": name,
        "mw": mw,
        "n_atoms": n_atoms,
        "n_bonds": n_bonds,
        "n_co": n_co,
        "co_ratio": n_co / n_bonds if n_bonds > 0 else 0,
        "rotatable": rot,
    }

    with open("experiments/iter_28_descriptors.json", "w") as f:
        json.dump([result], f, indent=2)

    print("\nDescriptors written to experiments/iter_28_descriptors.json")
else:
    print("Failed to parse SMILES")

print("\nDone!")
