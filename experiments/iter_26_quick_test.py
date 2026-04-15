"""Quick test of the SASA approach"""

import sys, os

print("STEP 1: Imports starting")
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

print("STEP 2: Loading PROTACs")
df = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t")
print(f"Loaded {len(df)} PROTACs")

print("\nSTEP 3: Processing first PROTAC")
row = df.iloc[0]
smiles = row["SMILES"]
name = row["Name"]
print(f"Name: {name}")
print(f"SMILES: {smiles[:50]}...")

print("\nSTEP 4: Parsing SMILES")
mol = Chem.MolFromSmiles(smiles)
print(f"Parsed molecule with {mol.GetNumAtoms()} heavy atoms")

print("\nSTEP 5: Adding hydrogens")
mol = Chem.AddHs(mol)
print(f"Now has {mol.GetNumAtoms()} total atoms")

print("\nSTEP 6: Generating conformers")
AllChem.EmbedMultipleConfs(mol, numConfs=5, randomSeed=42)
print(f"Generated {mol.GetNumConformers()} conformers")

print("\nSTEP 7: Minimizing")
for i in range(mol.GetNumConformers()):
    print(f"  Minimizing conf {i + 1}/{mol.GetNumConformers()}...")
    AllChem.MMFFOptimizeMolecule(mol, confId=i)
print("Minimization complete")

print("\nSTEP 8: Computing SASA")
from rdkit.Chem import rdFreeSASA

radii = rdFreeSASA.classifyAtoms(mol)
sasa_vals = []
for i in range(mol.GetNumConformers()):
    sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=i)
    sasa_vals.append(sasa)
    print(f"  Conf {i}: {sasa:.1f} A^2")

print(f"\nSASA mean: {np.mean(sasa_vals):.1f}")
print(f"SASA std: {np.std(sasa_vals):.1f}")

print("\nSUCCESS!")
