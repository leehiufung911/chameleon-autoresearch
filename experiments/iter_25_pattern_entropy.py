#!/usr/bin/env python3
"""
Iteration 25: Linker Pattern Entropy Analysis

Hypothesis: Chameleonic molecules have diverse linker patterns (mix of PEG and alkyl)
that create conformational heterogeneity, while non-chameleonic molecules have
homogeneous linkers (pure alkyl) that constrain to a single state.

Novel Approach: Compute entropy of linker composition patterns in the SMILES string.
High entropy = diverse linker = chameleonic. Low entropy = homogeneous = non-chameleonic.

This is a pure 2D topological approach that avoids the implicit solvent bias
by not generating any 3D conformers.
"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import json
import re
from collections import Counter
import math

output_lines = []


def log(msg):
    """Log to both stdout and collection."""
    output_lines.append(str(msg))
    print(str(msg), flush=True)


def compute_pattern_entropy(smiles, name):
    """
    Compute entropy of chemical patterns in the SMILES.

    The idea: chameleonic molecules have diverse patterns (mix of polar and
    non-polar segments), while non-chameleonic have homogeneous patterns.

    We extract n-grams from the SMILES and compute their entropy.
    """
    # Clean SMILES (remove stereo info for pattern matching)
    clean_smiles = re.sub(r"[@\[\]]", "", smiles)

    # Extract 3-character patterns (linker-relevant)
    patterns = []
    for i in range(len(clean_smiles) - 2):
        pattern = clean_smiles[i : i + 3]
        # Filter to only include patterns with C, O, N
        if any(c in pattern for c in "CON"):
            patterns.append(pattern)

    if not patterns:
        return None

    # Count pattern frequencies
    pattern_counts = Counter(patterns)
    total = len(patterns)

    # Compute Shannon entropy
    entropy = 0.0
    for count in pattern_counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)

    # Normalize by max possible entropy
    max_entropy = math.log2(len(pattern_counts)) if len(pattern_counts) > 1 else 1.0
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0.0

    # Count specific linker-related patterns
    peg_indicators = ["COC", "OCC", "CCO", "OCO"]
    alkyl_indicators = ["CCC", "CCCC"]
    amide_indicators = ["C(=O)N", "NC(=O)"]

    peg_count = sum(clean_smiles.count(p) for p in peg_indicators)
    alkyl_count = sum(clean_smiles.count(p) for p in alkyl_indicators)
    amide_count = sum(clean_smiles.count(p) for p in amide_indicators)

    # Compute pattern diversity ratio
    total_patterns = (
        peg_count + alkyl_count + amide_count + 1
    )  # +1 to avoid div by zero
    peg_ratio = peg_count / total_patterns
    alkyl_ratio = alkyl_count / total_patterns

    # Diversity score: high when both PEG and alkyl present
    diversity = 1 - abs(peg_ratio - alkyl_ratio)

    # Chameleonicity score
    chameleon_score = normalized_entropy * diversity * 10

    return {
        "name": name,
        "pattern_entropy": entropy,
        "normalized_entropy": normalized_entropy,
        "peg_count": peg_count,
        "alkyl_count": alkyl_count,
        "amide_count": amide_count,
        "peg_ratio": peg_ratio,
        "alkyl_ratio": alkyl_ratio,
        "diversity": diversity,
        "chameleon_score": chameleon_score,
    }


def main():
    log("=" * 80)
    log("ITERATION 25: Linker Pattern Entropy Analysis")
    log("=" * 80)
    log("")
    log("Hypothesis: Chameleonic molecules have diverse linker patterns")
    log("(high pattern entropy), while non-chameleonic have homogeneous patterns.")
    log("")
    log("Physical basis: PEG-alkyl mix creates conformational diversity")
    log("Pure alkyl constrains to extended state (low entropy)")
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
        result = compute_pattern_entropy(smiles, name)
        if result:
            results.append(result)
            log(f"\n{name}:")
            log(f"  Pattern entropy: {result['pattern_entropy']:.3f}")
            log(f"  Normalized entropy: {result['normalized_entropy']:.3f}")
            log(f"  PEG indicators: {result['peg_count']}")
            log(f"  Alkyl indicators: {result['alkyl_count']}")
            log(f"  Diversity: {result['diversity']:.3f}")
            log(f"  Chameleon score: {result['chameleon_score']:.3f}")

    # Separation check
    log("\n" + "-" * 80)
    log("SEPARATION CHECK")
    log("-" * 80)

    cham_scores = [
        r["chameleon_score"] for r in results if r["name"] in ["protac_1", "protac_2"]
    ]
    non_scores = [r["chameleon_score"] for r in results if r["name"] == "protac_3"]

    if cham_scores and non_scores:
        cham_mean = sum(cham_scores) / len(cham_scores)
        non_mean = sum(non_scores) / len(non_scores)
        log(f"\nChameleonic avg: {cham_mean:.3f}")
        log(f"Non-chameleonic: {non_mean:.3f}")
        log(f"Difference: {cham_mean - non_mean:.3f}")
        if cham_mean > non_mean:
            log("Status: CORRECT ORDER (chameleonic > non-chameleonic)")
        else:
            log("Status: REVERSED ORDER")

    # Load benchmark - try both files
    log("\n" + "-" * 80)
    log("BENCHMARK DESCRIPTORS")
    log("-" * 80)

    benchmark = []

    # Try labelled_set.tsv first
    try:
        with open("chameleon_local/labelled_set.tsv", "r") as f:
            lines = f.readlines()
            log(f"Read {len(lines)} lines from labelled_set.tsv")
            for i, line in enumerate(lines[1:]):  # Skip header
                parts = line.strip().split("\t")
                log(f"Line {i + 1}: {len(parts)} parts")
                if len(parts) >= 3:
                    benchmark.append((parts[0], parts[1], parts[2]))
    except Exception as e:
        log(f"Error reading labelled_set.tsv: {e}")

    if not benchmark:
        # Fall back to benchmark.tsv
        try:
            with open("chameleon_local/benchmark.tsv", "r") as f:
                lines = f.readlines()
                log(f"Read {len(lines)} lines from benchmark.tsv")
                for i, line in enumerate(lines):
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        # No label in this file, default to 'unknown'
                        benchmark.append((parts[0], parts[1], "unknown"))
        except Exception as e:
            log(f"Error reading benchmark.tsv: {e}")

    log(f"\nLoaded {len(benchmark)} molecules from benchmark")

    # Process all benchmark molecules
    benchmark_results = []
    errors = []
    for i, (name, smiles, label) in enumerate(benchmark):
        try:
            result = compute_pattern_entropy(smiles, name)
            if result:
                result["label"] = label
                benchmark_results.append(result)
        except Exception as e:
            errors.append((name, str(e)))

    log(f"\nProcessed {len(benchmark_results)} benchmark molecules")
    if errors:
        log(f"Errors: {len(errors)}")
        for name, err in errors[:3]:
            log(f"  {name}: {err}")

    # Statistics
    cham_results = [r for r in benchmark_results if r["label"] == "cham"]
    non_results = [r for r in benchmark_results if r["label"] == "non"]

    log(f"\nChameleonic: {len(cham_results)} molecules")
    log(f"Non-chameleonic: {len(non_results)} molecules")

    if cham_results and non_results:
        cham_scores = [r["chameleon_score"] for r in cham_results]
        non_scores = [r["chameleon_score"] for r in non_results]

        import numpy as np

        log(
            f"\nChameleonic mean: {np.mean(cham_scores):.3f} +/- {np.std(cham_scores):.3f}"
        )
        log(
            f"Non-chameleonic mean: {np.mean(non_scores):.3f} +/- {np.std(non_scores):.3f}"
        )
        log(f"Difference: {np.mean(cham_scores) - np.mean(non_scores):.3f}")

        # Top and bottom scorers
        all_scores = [
            (r["name"], r["chameleon_score"], r["label"]) for r in benchmark_results
        ]
        all_scores.sort(key=lambda x: x[1], reverse=True)

        log(f"\nTop 5 chameleonic scores:")
        for name, score, label in all_scores[:5]:
            log(f"  {name}: {score:.3f} ({label})")

        log(f"\nBottom 5 chameleonic scores:")
        for name, score, label in all_scores[-5:]:
            log(f"  {name}: {score:.3f} ({label})")

    log("\n" + "=" * 80)
    log("END OF ANALYSIS")
    log("=" * 80)

    # Write output
    output_text = "\n".join(output_lines)
    with open("experiments/iter_25_output.txt", "w") as f:
        f.write(output_text)

    # Save descriptors - include all with proper structure
    descriptors = []
    for r in benchmark_results:
        descriptors.append(
            {
                "name": r["name"],
                "label": r["label"],
                "chameleon_score": r["chameleon_score"],
                "pattern_entropy": r["pattern_entropy"],
                "peg_count": r["peg_count"],
                "alkyl_count": r["alkyl_count"],
            }
        )

    with open("experiments/iter_25_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)

    log(f"\nResults written to experiments/iter_25_output.txt")
    log(f"Descriptors written to experiments/iter_25_descriptors.json")


if __name__ == "__main__":
    main()
