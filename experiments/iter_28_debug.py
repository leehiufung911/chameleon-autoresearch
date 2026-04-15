import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

print("DEBUG: Starting script")
sys.stderr.write("DEBUG: stderr write\n")

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

print("DEBUG: Imports successful")

# Test loading PROTACs
with open("chameleon_local/user_protacs.tsv", "r") as f:
    lines = f.readlines()
    print(f"DEBUG: Loaded {len(lines)} PROTACs")
    for line in lines:
        parts = line.strip().split(maxsplit=1)
        print(f"DEBUG: {parts[0]}")

# Try to process one
smiles = lines[0].strip().split(maxsplit=1)[1]
print(f"DEBUG: SMILES length {len(smiles)}")

mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("DEBUG: Failed to parse SMILES")
else:
    print(f"DEBUG: Parsed molecule with {mol.GetNumAtoms()} atoms")

mol = Chem.AddHs(mol)
print(f"DEBUG: Added Hs, now {mol.GetNumAtoms()} atoms")

# Try conformer generation
ps = AllChem.ETKDGv3()
ps.randomSeed = 42
AllChem.EmbedMultipleConfs(mol, 5, ps)
print(f"DEBUG: Generated {mol.GetNumConformers()} conformers")

print("DEBUG: Script completed successfully")
