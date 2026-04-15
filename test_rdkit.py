#!/usr/bin/env python3
import sys

print("Step 1: Starting...", flush=True)
import numpy as np

print("Step 2: NumPy imported", flush=True)
from rdkit import Chem

print("Step 3: RDKit imported", flush=True)
mol = Chem.MolFromSmiles("CCC")
print("Step 4: Molecule created:", mol is not None, flush=True)
print("DONE", flush=True)
