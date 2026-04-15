print("START")
import sys

sys.path.insert(0, "chameleon_local")
print("PATH SET")

from rdkit import Chem

print("RDKIT IMPORTED")

with open("chameleon_local/user_protacs.tsv") as f:
    lines = f.readlines()
print(f"LINES: {len(lines)}")

parts = lines[0].strip().split("\t")
print(f"NAME: {parts[0]}")

mol = Chem.MolFromSmiles(parts[1])
print(f"ATOMS: {mol.GetNumAtoms()}")

with open("experiments/iter_19_output.txt", "w") as f:
    f.write(f"Successfully parsed {len(lines)} PROTACs\n")
    f.write(f"First molecule: {parts[0]} with {mol.GetNumAtoms()} atoms\n")

print("DONE")
