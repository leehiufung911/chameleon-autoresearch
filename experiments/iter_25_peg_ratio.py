#!/usr/bin/env python3
"""
Iteration 25: Linker PEG Ratio Analysis

Hypothesis: Chameleonicity is driven by PEG linker content, not pattern diversity.
PEG linkers (gauche preference) enable folding; alkyl linkers (anti preference) resist it.

Novel Approach: Compute PEG/(PEG+alkyl) ratio. High ratio = chameleonic.
Pure 2D approach - no 3D conformers.
"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import json
import numpy as np

output_lines = []


def log(msg):
    output_lines.append(str(msg))
    print(str(msg), flush=True)


def compute_peg_ratio(smiles, name):
    """Compute PEG ratio from SMILES patterns."""
    # Count PEG patterns (C-O-C)
    peg = smiles.count("COC") + smiles.count("OCCO") + smiles.count("CCOCC")

    # Count alkyl patterns (C-C-C)
    alkyl = smiles.count("CCC") + smiles.count("CCCC") + smiles.count("CCCCC")

    # Avoid division by zero
    total = peg + alkyl + 1

    # PEG ratio
    peg_ratio = peg / total
    alkyl_ratio = alkyl / total

    # Chameleonicity score
    # High peg_ratio = chameleonic
    # Normalize by molecule size
    score = peg_ratio * 10

    return {
        "name": name,
        "peg_count": peg,
        "alkyl_count": alkyl,
        "peg_ratio": peg_ratio,
        "alkyl_ratio": alkyl_ratio,
        "chameleon_score": score,
    }


def main():
    log("=" * 80)
    log("ITERATION 25: PEG Linker Ratio Analysis")
    log("=" * 80)
    log("")
    log("Hypothesis: Chameleonicity correlates with PEG content in linker.")
    log("PEG (gauche) enables folding; alkyl (anti) resists it.")
    log("")

    # Load PROTACs
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
        result = compute_peg_ratio(smiles, name)
        if result:
            results.append(result)
            log(f"\n{name}:")
            log(f"  PEG patterns: {result['peg_count']}")
            log(f"  Alkyl patterns: {result['alkyl_count']}")
            log(f"  PEG ratio: {result['peg_ratio']:.3f}")
            log(f"  Chameleon score: {result['chameleon_score']:.3f}")

    # Separation check
    log("\n" + "-" * 80)
    log("SEPARATION CHECK")
    log("-" * 80)

    cham = [r for r in results if r["name"] in ["protac_1", "protac_2"]]
    non = [r for r in results if r["name"] == "protac_3"]

    if cham and non:
        cham_mean = np.mean([r["chameleon_score"] for r in cham])
        non_mean = np.mean([r["chameleon_score"] for r in non])
        log(f"\nChameleonic avg: {cham_mean:.3f}")
        log(f"Non-chameleonic: {non_mean:.3f}")
        log(f"Difference: {cham_mean - non_mean:.3f}")
        if cham_mean > non_mean:
            log("Status: CORRECT ORDER")
        else:
            log("Status: REVERSED ORDER")

    # Load benchmark
    log("\n" + "-" * 80)
    log("BENCHMARK DESCRIPTORS")
    log("-" * 80)

    benchmark = []
    with open("chameleon_local/labelled_set.tsv", "r") as f:
        for line in f.readlines()[1:]:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                benchmark.append((parts[0], parts[1], parts[2]))

    log(f"Loaded {len(benchmark)} molecules")

    # Process benchmark
    bench_results = []
    for name, smiles, label in benchmark:
        result = compute_peg_ratio(smiles, name)
        if result:
            result["label"] = label
            bench_results.append(result)

    log(f"\nProcessed {len(bench_results)} molecules")

    # Statistics
    cham_results = [r for r in bench_results if r["label"] == "cham"]
    non_results = [r for r in bench_results if r["label"] == "non"]

    log(f"\nChameleonic: {len(cham_results)}")
    log(f"Non-chameleonic: {len(non_results)}")

    if cham_results and non_results:
        cham_scores = [r["chameleon_score"] for r in cham_results]
        non_scores = [r["chameleon_score"] for r in non_results]

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

    # Save descriptors
    with open("experiments/iter_25_descriptors.json", "w") as f:
        json.dump(bench_results, f, indent=2)

    log(f"\nWritten to experiments/iter_25_output.txt and iter_25_descriptors.json")


if __name__ == "__main__":
    main()
