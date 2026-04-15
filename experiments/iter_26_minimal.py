"""Minimal SASA test"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"

print("Starting minimal test...", flush=True)

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA

print("RDKit imports successful", flush=True)

# Simple test molecule
mol = Chem.MolFromSmiles("CCOCC")
mol = Chem.AddHs(mol)
print(f"Molecule: {mol.GetNumAtoms()} atoms", flush=True)

# Embed
AllChem.EmbedMultipleConfs(mol, numConfs=5, randomSeed=42)
print(f"Embedded: {mol.GetNumConformers()} conformers", flush=True)

# Minimize
for i in range(mol.GetNumConformers()):
    AllChem.MMFFOptimizeMolecule(mol, confId=i)
print("Minimization complete", flush=True)

# SASA
radii = rdFreeSASA.classifyAtoms(mol)
sasa_vals = []
for i in range(mol.GetNumConformers()):
    sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=i)
    sasa_vals.append(sasa)
    print(f"  Conf {i}: SASA = {sasa:.1f} A^2", flush=True)

print(f"\nSASA mean: {np.mean(sasa_vals):.1f}", flush=True)
print(f"SASA std: {np.std(sasa_vals):.1f}", flush=True)
print("SUCCESS!", flush=True)
