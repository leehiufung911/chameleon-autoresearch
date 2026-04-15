import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

print("=" * 70)
print("ITERATION 28: SMILES-Based Topology Analysis")
print("=" * 70)
print()

from rdkit import Chem
from rdkit.Chem import Descriptors
import json

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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue

    mol = Chem.AddHs(mol)
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()
    mw = Descriptors.MolWt(mol)
    rot = Descriptors.NumRotatableBonds(mol)

    # Count C-O bonds (simple ethers)
    n_co = 0
    for bond in mol.GetBonds():
        n1 = bond.GetBeginAtom().GetAtomicNum()
        n2 = bond.GetEndAtom().GetAtomicNum()
        if (n1 == 6 and n2 == 8) or (n1 == 8 and n2 == 6):
            n_co += 1

    label = labels.get(name, "?")
    results.append(
        {
            "name": name,
            "label": label,
            "mw": mw,
            "n_atoms": n_atoms,
            "n_bonds": n_bonds,
            "n_co": n_co,
            "co_ratio": n_co / n_bonds if n_bonds > 0 else 0,
            "rotatable": rot,
        }
    )

    print(
        f"{name} ({label}): MW={mw:.1f}, CO_bonds={n_co}, CO_ratio={n_co / n_bonds:.3f}, rot={rot}"
    )

# Summary
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)

cham_co = [r["co_ratio"] for r in results if r["label"] == "CHAM"]
non_co = [r["co_ratio"] for r in results if r["label"] == "NON"]

if cham_co and non_co:
    cham_mean = sum(cham_co) / len(cham_co)
    non_mean = sum(non_co) / len(non_co)
    gap = cham_mean - non_mean

    print(f"Chameleonic avg CO ratio: {cham_mean:.4f}")
    print(f"Non-chameleonic avg: {non_mean:.4f}")
    print(f"Separation: {gap:+.4f}")

    if gap > 0.01:
        print("Verdict: SUCCESS")
    elif gap > 0:
        print("Verdict: PARTIAL SUCCESS")
    else:
        print("Verdict: FAILED")

# Write descriptors
with open("experiments/iter_28_descriptors.json", "w") as f:
    json.dump(results, f, indent=2)

print()
print("Descriptors written to experiments/iter_28_descriptors.json")
