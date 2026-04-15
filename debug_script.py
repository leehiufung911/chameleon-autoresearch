#!/usr/bin/env python3
import sys

sys.stderr.write("DEBUG: Script started\n")
sys.stderr.flush()

import os

os.environ["PYTHONUNBUFFERED"] = "1"

sys.stderr.write("DEBUG: Importing numpy...\n")
sys.stderr.flush()
import numpy as np

sys.stderr.write("DEBUG: Importing RDKit...\n")
sys.stderr.flush()
from rdkit import Chem

sys.stderr.write("DEBUG: Importing Descriptors...\n")
sys.stderr.flush()
from rdkit.Chem import Descriptors

sys.stderr.write("DEBUG: Testing SMILES...\n")
sys.stderr.flush()
mol = Chem.MolFromSmiles("CCC")
sys.stderr.write(f"DEBUG: Molecule created: {mol is not None}\n")
sys.stderr.flush()

sys.stderr.write("DEBUG: Script complete!\n")
sys.stderr.flush()
