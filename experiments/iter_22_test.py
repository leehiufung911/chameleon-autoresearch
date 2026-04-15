import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys

sys.path.insert(0, "chameleon_local")

print("Step 1: Importing RDKit...")
from rdkit import Chem

print("Step 2: RDKit imported successfully")

from rdkit.Chem import rdMolDescriptors

print("Step 3: rdMolDescriptors imported")

print("Step 4: Loading user_protacs.tsv...")
protacs = []
with open("chameleon_local/user_protacs.tsv", "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split(" ", 1)
        if len(parts) == 2:
            name = parts[0]
            rest = parts[1]
            tokens = rest.split()
            smiles = tokens[0] if tokens else rest
            label = "CHAM" if name in ["protac_1", "protac_2"] else "NON"
            print(f"  Loading {name}: {smiles[:30]}...")
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                protacs.append(
                    {"name": name, "mol": mol, "smiles": smiles, "label": label}
                )
                print(f"    Success: {mol.GetNumAtoms()} atoms")
            else:
                print(f"    FAILED to parse SMILES")

print(f"\nStep 5: Loaded {len(protacs)} PROTACs")

print("\nStep 6: Testing fragment counting...")
for item in protacs:
    mol = item["mol"]
    name = item["name"]

    # PEG pattern
    peg_pattern = Chem.MolFromSmarts("[OX2][C;!R][C;!R][OX2]")
    peg_matches = mol.GetSubstructMatches(peg_pattern)

    # Amide pattern
    amide_pattern = Chem.MolFromSmarts("[NX3](C=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    # Long alkyl
    alkyl_pattern = Chem.MolFromSmarts("[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]")
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)

    print(
        f"{name}: PEG={len(peg_matches)}, Amide={len(amide_matches)}, Alkyl5={len(alkyl_matches)}"
    )

print("\nStep 7: Writing output...")
output = "Test completed successfully!\n"
with open("experiments/iter_22_output.txt", "w") as f:
    f.write(output)
print(output)
