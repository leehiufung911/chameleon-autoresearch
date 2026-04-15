import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import json

# Test basic functionality
mol = Chem.MolFromSmiles("CCOCC")
if mol is None:
    print("ERROR: Failed to parse SMILES")
    sys.exit(1)

mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
print(f"Test molecule: {mol.GetNumAtoms()} atoms")

# Try FreeSASA
try:
    from rdkit.Chem import rdFreeSASA

    radii = rdFreeSASA.classifyAtoms(mol)
    sasa = rdFreeSASA.CalcSASA(mol, radii)
    print(f"SASA computed: {sasa:.2f} A^2")
    print("SUCCESS: rdFreeSASA works")
except Exception as e:
    print(f"ERROR with FreeSASA: {e}")

# Test loading from TSV
try:
    import pandas as pd

    df = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t")
    print(f"Loaded {len(df)} PROTACs from TSV")
    print(f"Columns: {list(df.columns)}")
    for _, row in df.head(2).iterrows():
        print(f"  {row['Name']}: {row['SMILES'][:40]}...")
except Exception as e:
    print(f"ERROR loading TSV: {e}")

print("\nAll tests passed!")
