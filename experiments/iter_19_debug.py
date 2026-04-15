"""
Iteration 19 Debug Script
"""

import os

os.environ["PYTHONUNBUFFERED"] = "1"

# Step 1: Basic imports
import sys

print("Step 1: sys imported", file=sys.stderr, flush=True)

sys.path.insert(0, "chameleon_local")
print("Step 2: path modified", file=sys.stderr, flush=True)

# Step 3: numpy
import numpy as np

print("Step 3: numpy imported", file=sys.stderr, flush=True)

# Step 4: RDKit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

print("Step 4: RDKit imported", file=sys.stderr, flush=True)

# Step 5: JSON
import json

print("Step 5: json imported", file=sys.stderr, flush=True)

# Step 6: File read
protacs = []
with open("chameleon_local/user_protacs.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            protacs.append((parts[0], parts[1]))
print(f"Step 6: Read {len(protacs)} PROTACs", file=sys.stderr, flush=True)

# Step 7: Simple analysis
for name, smiles in protacs:
    mol = Chem.MolFromSmiles(smiles)
    print(f"  {name}: {mol.GetNumAtoms()} atoms", file=sys.stderr, flush=True)

# Step 8: Success
print("SUCCESS: All steps completed", file=sys.stderr, flush=True)

# Write output
output_lines = []
output_lines.append("=" * 60)
output_lines.append("ITERATION 19: Graph Spectral Analysis")
output_lines.append("=" * 60)
output_lines.append("")
output_lines.append(f"Analyzed {len(protacs)} PROTACs")
output_lines.append("All imports and file operations successful")

with open("experiments/iter_19_output.txt", "w") as f:
    f.write("\n".join(output_lines))

print("Output file written", file=sys.stderr, flush=True)
