import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"

print("Test 1: Starting script", file=sys.stderr)

from rdkit import Chem

print("Test 2: RDKit imported", file=sys.stderr)

with open("chameleon_local/user_protacs.tsv", "r") as f:
    lines = f.readlines()
    print(f"Test 3: Read {len(lines)} lines", file=sys.stderr)

for i, line in enumerate(lines):
    parts = line.strip().split()
    print(f"Test 4: Line {i} parts: {len(parts)}", file=sys.stderr)

print("Test 5: Done", file=sys.stderr)
