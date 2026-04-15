import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"

log = open("experiments/iter_25_rdkit_test.log", "w")
log.write("Step 1: Start\n")
log.flush()

log.write("Step 2: Importing Chem...\n")
log.flush()
from rdkit import Chem

log.write("Step 3: Chem imported\n")
log.flush()

log.write("Step 4: Importing AllChem...\n")
log.flush()
from rdkit.Chem import AllChem

log.write("Step 5: AllChem imported\n")
log.flush()

log.write("Step 6: Parsing simple SMILES...\n")
log.flush()
mol = Chem.MolFromSmiles("CCO")
log.write(f"Step 7: Parsed molecule with {mol.GetNumAtoms()} atoms\n")
log.flush()

log.write("Step 8: Adding Hs...\n")
log.flush()
mol = Chem.AddHs(mol)
log.write(f"Step 9: Now has {mol.GetNumAtoms()} atoms\n")
log.flush()

log.write("Step 10: Generating conformers...\n")
log.flush()
try:
    AllChem.EmbedMultipleConfs(mol, 3, AllChem.ETKDGv3())
    log.write(f"Step 11: Generated {mol.GetNumConformers()} conformers\n")
    log.flush()
except Exception as e:
    log.write(f"Step 11 FAILED: {e}\n")
    log.flush()

log.write("Step 12: Done\n")
log.flush()
log.close()

print("Script completed", file=sys.stderr)
