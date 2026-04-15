import os

os.environ["PYTHONUNBUFFERED"] = "1"

print("Starting minimal test")

from rdkit import Chem

print("Imported Chem")

# Don't import AllChem separately - use Chem.AllChem
mol = Chem.MolFromSmiles("CCCC")
print("Created mol")

mol = Chem.AddHs(mol)
print("Added Hs")

params = Chem.AllChem.ETKDGv3()
print("Created params")

cids = Chem.AllChem.EmbedMultipleConfs(mol, numConfs=2, params=params)
print(f"Generated {len(cids)} conformers")

if len(cids) > 0:
    Chem.AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=10, mmffVariant="MMFF94s")
    print("Optimized")

print("Success!")
