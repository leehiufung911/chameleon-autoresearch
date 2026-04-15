import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys

sys.path.insert(0, "chameleon_local")

log = []
log.append("=== Iteration 19 Spectral Analysis ===")
log.append("")

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import json

log.append("Imports successful")


def build_weighted_graph(mol):
    n_atoms = mol.GetNumAtoms()
    adj = np.zeros((n_atoms, n_atoms))
    bond_weight = {
        Chem.BondType.SINGLE: 1.0,
        Chem.BondType.DOUBLE: 2.0,
        Chem.BondType.TRIPLE: 3.0,
        Chem.BondType.AROMATIC: 1.5,
    }

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        wt = bond_weight.get(bond.GetBondType(), 1.0)
        atom_i = mol.GetAtom(i)
        atom_j = mol.GetAtom(j)
        polar_boost = (
            1.3
            if atom_i.GetSymbol() in ["N", "O"] or atom_j.GetSymbol() in ["N", "O"]
            else 1.0
        )
        adj[i, j] = wt * polar_boost
        adj[j, i] = wt * polar_boost

    return adj


def compute_laplacian(adj):
    degrees = np.sum(adj, axis=1)
    degrees = np.where(degrees == 0, 1, degrees)
    d_inv_sqrt = np.diag(1.0 / np.sqrt(degrees))
    return np.eye(len(adj)) - d_inv_sqrt @ adj @ d_inv_sqrt


def compute_spectral(mol):
    adj = build_weighted_graph(mol)
    laplacian = compute_laplacian(adj)

    try:
        eigenvalues = np.linalg.eigvalsh(laplacian)
        eigenvalues = np.sort(eigenvalues)
    except:
        return None

    non_zero = eigenvalues[eigenvalues > 1e-10]
    if len(non_zero) < 2:
        return None

    fiedler = non_zero[0]
    eig_probs = non_zero / np.sum(non_zero)
    entropy = -np.sum(eig_probs * np.log(eig_probs + 1e-10))
    resistance = np.sum(1.0 / non_zero) / (mol.GetNumAtoms() ** 2)

    return {
        "fiedler": fiedler,
        "spectral_entropy": entropy,
        "normalized_resistance": resistance,
    }


def analyze_mol(smiles, name):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    desc = compute_spectral(mol)
    if desc is None:
        return None

    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Compute diffusion score
    if n_rot > 0:
        flex = np.log1p(n_rot)
        score = (
            desc["fiedler"] * 2.0
            + (5.0 - desc["spectral_entropy"]) * 0.5
            + (1.0 / (desc["normalized_resistance"] + 0.01)) * 0.3
        ) / flex
    else:
        score = desc["fiedler"] * 2.0

    return {
        "name": name,
        "smiles": smiles,
        "diffusion_score": score,
        "n_rot": n_rot,
        **desc,
    }


# Read PROTACs
log.append("Reading PROTACs...")
protacs = []
with open("chameleon_local/user_protacs.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            protacs.append((parts[0], parts[1]))

log.append(f"Loaded {len(protacs)} PROTACs")
log.append("")

# Analyze PROTACs
results = []
for name, smiles in protacs:
    res = analyze_mol(smiles, name)
    if res:
        results.append(res)
        gt = "CHAM" if name in ["protac_1", "protac_2"] else "NON"
        log.append(f"{name} ({gt}):")
        log.append(f"  Fiedler: {res['fiedler']:.4f}")
        log.append(f"  Entropy: {res['spectral_entropy']:.3f}")
        log.append(f"  Resistance: {res['normalized_resistance']:.4f}")
        log.append(f"  Score: {res['diffusion_score']:.4f}")
        log.append("")

# Check separation
cham_scores = [
    r["diffusion_score"] for r in results if r["name"] in ["protac_1", "protac_2"]
]
non_scores = [r["diffusion_score"] for r in results if r["name"] == "protac_3"]

if cham_scores and non_scores:
    log.append("SEPARATION ANALYSIS")
    log.append(f"Chameleonic avg: {np.mean(cham_scores):.4f}")
    log.append(f"Non-chameleonic avg: {np.mean(non_scores):.4f}")
    log.append(f"Difference: {np.mean(cham_scores) - np.mean(non_scores):+.4f}")
    log.append("")

    if np.mean(cham_scores) > np.mean(non_scores):
        log.append("RESULT: CORRECT ORDERING - SUCCESS")
    else:
        log.append("RESULT: FAILED - Wrong ordering")
    log.append("")

# Analyze benchmark
log.append("BENCHMARK ANALYSIS (49 molecules)")
benchmark = []
with open("chameleon_local/benchmark.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            benchmark.append((parts[0], parts[1]))

log.append(f"Loaded {len(benchmark)} benchmark molecules")
log.append("")

benchmark_results = []
for name, smiles in benchmark:
    res = analyze_mol(smiles, name)
    if res:
        benchmark_results.append(res)

log.append(f"Successfully analyzed {len(benchmark_results)} molecules")
log.append("")

# Summary
log.append("=" * 60)
log.append("SUMMARY")
log.append("=" * 60)
log.append("Graph Spectral Diffusion Analysis from network science:")
log.append("- Fiedler value: Graph connectivity measure")
log.append("- Spectral entropy: Folding landscape order")
log.append("- Effective resistance: Pathway efficiency")
log.append("")
log.append("This is a purely 2D topological approach.")
log.append("No 3D conformers, no energy calculations.")

# Write output
output = "\n".join(log)
with open("experiments/iter_19_output.txt", "w") as f:
    f.write(output)

with open("experiments/iter_19_descriptors.json", "w") as f:
    json.dump(benchmark_results, f, indent=2)

print("SUCCESS: Wrote output files", flush=True)
