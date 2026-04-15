import sys, os
os.environ['PYTHONUNBUFFERED'] = '1'
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem

print("Testing imports...")
print(f"NumPy version: {np.__version__}")
print(f"RDKit available: {Chem is not None}")

# Test SMILES parsing
smiles = "CCO"
mol = Chem.MolFromSmiles(smiles)
print(f"Parsed {smiles}: {mol is not None}")

print("SUCCESS: All imports working")
