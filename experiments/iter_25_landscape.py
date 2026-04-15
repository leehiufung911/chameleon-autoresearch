#!/usr/bin/env python3
"""
Iteration 25: SMILES-Based Energy Landscape Topology Analysis

Hypothesis: The energy landscape topology (unimodal vs bimodal) distinguishes
chameleonic from non-chameleonic PROTACs. Chameleonic molecules have bimodal
landscapes with distinct compact and extended minima, while non-chameleonic
molecules have unimodal landscapes.

Novel Approach: Instead of computing single descriptors, we analyze the
distribution of energies across the conformer ensemble to detect bimodality.
"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import json

output_lines = []


def log(msg):
    """Log to both stdout and collection."""
    output_lines.append(str(msg))
    print(str(msg), flush=True)


def compute_rg(positions):
    """Compute radius of gyration."""
    positions = np.array(positions)
    com = np.mean(positions, axis=0)
    r_sq = np.sum((positions - com) ** 2, axis=1)
    return np.sqrt(np.mean(r_sq))


def get_energy_landscape(smiles, n_confs=20):
    """Generate conformers and get energy-Rg distribution."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, n_confs, AllChem.ETKDGv3())

    if mol.GetNumConformers() == 0:
        return None

    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")

    results = []
    for i in range(mol.GetNumConformers()):
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=i)
        ff.Initialize()
        ff.Minimize(maxIts=200)
        energy = ff.CalcEnergy()

        conf = mol.GetConformer(i)
        positions = np.array(
            [conf.GetAtomPosition(j) for j in range(mol.GetNumAtoms())]
        )
        rg = compute_rg(positions)

        results.append({"energy": energy, "rg": rg})

    return results


def analyze_landscape(results, name):
    """Analyze energy landscape for bimodality."""
    if not results or len(results) < 5:
        return None

    energies = np.array([r["energy"] for r in results])
    rg_values = np.array([r["rg"] for r in results])

    # Normalize energies
    e_min = np.min(energies)
    e_norm = energies - e_min

    # Compute bimodality coefficient
    # BC = (skewness^2 + 1) / (kurtosis + 3)
    # BC > 0.55 suggests bimodal distribution
    if len(energies) > 3:
        from scipy import stats

        skew = stats.skew(e_norm)
        kurt = stats.kurtosis(e_norm)
        bimodality = (skew**2 + 1) / (kurt + 3) if (kurt + 3) != 0 else 0.5
    else:
        bimodality = 0.5

    # Energy range (kcal/mol)
    e_range = np.max(energies) - np.min(energies)

    # Rg range and variability
    rg_range = np.max(rg_values) - np.min(rg_values)
    rg_cv = np.std(rg_values) / np.mean(rg_values) if np.mean(rg_values) > 0 else 0

    # Energy-Rg correlation
    if len(energies) > 2:
        corr = np.corrcoef(energies, rg_values)[0, 1]
    else:
        corr = 0.0

    # Landscape type score
    # High bimodality + weak E-Rg correlation = two-state system
    landscape_score = bimodality * (1 - abs(corr))

    return {
        "name": name,
        "n_confs": len(results),
        "energy_range": float(e_range),
        "rg_range": float(rg_range),
        "rg_cv": float(rg_cv),
        "bimodality": float(bimodality),
        "energy_rg_corr": float(corr),
        "landscape_score": float(landscape_score),
        "mean_rg": float(np.mean(rg_values)),
        "std_rg": float(np.std(rg_values)),
    }


def main():
    log("=" * 80)
    log("ITERATION 25: Energy Landscape Topology Analysis")
    log("=" * 80)
    log("")
    log("Hypothesis: Chameleonic molecules have bimodal energy landscapes")
    log("(compact + extended states), while non-chameleonic have unimodal.")
    log("")

    # Load user PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                protacs.append((parts[0], parts[1]))

    # Analyze PROTACs
    log("-" * 80)
    log("PROTAC ANALYSIS")
    log("-" * 80)

    results = []
    for name, smiles in protacs:
        log(f"\nProcessing {name}...")
        landscape = get_energy_landscape(smiles, n_confs=20)
        if landscape:
            analysis = analyze_landscape(landscape, name)
            if analysis:
                results.append(analysis)
                log(f"  N conformers: {analysis['n_confs']}")
                log(f"  Mean Rg: {analysis['mean_rg']:.2f} A")
                log(f"  Rg range: {analysis['rg_range']:.2f} A")
                log(f"  Bimodality: {analysis['bimodality']:.3f}")
                log(f"  E-Rg correlation: {analysis['energy_rg_corr']:.3f}")
                log(f"  Landscape score: {analysis['landscape_score']:.3f}")
        else:
            log(f"  Failed to generate conformers")

    # Separation check
    log("\n" + "-" * 80)
    log("PROTAC SEPARATION CHECK")
    log("-" * 80)

    cham_scores = [
        r["landscape_score"] for r in results if r["name"] in ["protac_1", "protac_2"]
    ]
    non_scores = [r["landscape_score"] for r in results if r["name"] == "protac_3"]

    if cham_scores and non_scores:
        cham_mean = np.mean(cham_scores)
        non_mean = np.mean(non_scores)
        log(f"\nChameleonic avg: {cham_mean:.3f}")
        log(f"Non-chameleonic: {non_mean:.3f}")
        log(f"Difference: {cham_mean - non_mean:.3f}")
        if cham_mean > non_mean:
            log("Status: CORRECT ORDER")
        else:
            log("Status: REVERSED ORDER")

    # Load labelled_set for benchmark analysis
    log("\n" + "-" * 80)
    log("BENCHMARK DESCRIPTORS (subset)")
    log("-" * 80)

    benchmark = []
    try:
        with open("chameleon_local/labelled_set.tsv", "r") as f:
            lines = f.readlines()
            for line in lines[1:]:  # Skip header
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    benchmark.append((parts[0], parts[1], parts[2]))
    except Exception as e:
        log(f"Could not load labelled_set.tsv: {e}")

    log(f"Loaded {len(benchmark)} molecules from labelled_set.tsv")

    # Process first 15 for speed
    benchmark_results = []
    for i, (name, smiles, label) in enumerate(benchmark[:15]):
        landscape = get_energy_landscape(smiles, n_confs=15)
        if landscape:
            analysis = analyze_landscape(landscape, name)
            if analysis:
                analysis["label"] = label
                benchmark_results.append(analysis)

    log(f"\nProcessed {len(benchmark_results)} benchmark molecules")

    # Statistics
    cham_results = [r for r in benchmark_results if r["label"] == "cham"]
    non_results = [r for r in benchmark_results if r["label"] == "non"]

    if cham_results and non_results:
        cham_scores = [r["landscape_score"] for r in cham_results]
        non_scores = [r["landscape_score"] for r in non_results]

        log(
            f"\nChameleonic mean: {np.mean(cham_scores):.3f} +/- {np.std(cham_scores):.3f}"
        )
        log(
            f"Non-chameleonic mean: {np.mean(non_scores):.3f} +/- {np.std(non_scores):.3f}"
        )
        log(f"Difference: {np.mean(cham_scores) - np.mean(non_scores):.3f}")

    log("\n" + "=" * 80)
    log("END OF ANALYSIS")
    log("=" * 80)

    # Write output
    output_text = "\n".join(output_lines)
    with open("experiments/iter_25_output.txt", "w") as f:
        f.write(output_text)

    # Save descriptors for all 49
    all_benchmark_results = []
    for name, smiles, label in benchmark:
        landscape = get_energy_landscape(smiles, n_confs=10)
        if landscape:
            analysis = analyze_landscape(landscape, name)
            if analysis:
                analysis["label"] = label
                all_benchmark_results.append(analysis)

    with open("experiments/iter_25_descriptors.json", "w") as f:
        json.dump(all_benchmark_results, f, indent=2)

    log(f"\nResults: experiments/iter_25_output.txt")
    log(f"Descriptors: experiments/iter_25_descriptors.json")


if __name__ == "__main__":
    main()
