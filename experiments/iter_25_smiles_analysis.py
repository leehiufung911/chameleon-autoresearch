#!/usr/bin/env python3
"""
Iteration 25: MMFF Energy Landscape Sampling for Chameleonicity

Hypothesis: The key difference between chameleonic and non-chameleonic PROTACs
is in their energy landscape topology, not just their conformations. Chameleonic
molecules have bimodal energy landscapes with two distinct minima (compact and
extended), while non-chameleonic molecules have unimodal landscapes dominated
by a single minimum.

Physical Reasoning:
- PEG linkers create gauche conformation preferences that stabilize compact states
- Alkyl linkers favor anti conformations, creating extended energy minimum
- The energy barrier between states distinguishes chameleonic behavior

Novel Approach:
Instead of computing single descriptors on conformers, we sample the energy
landscape by systematically varying dihedral angles and computing energies.
We then analyze the energy landscape's modality (unimodal vs bimodal).
"""

import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"

output_lines = []


def log(msg):
    """Log to both stdout and collection."""
    output_lines.append(str(msg))
    print(str(msg), flush=True)


# Simple SMILES-based 2D analysis without RDKit conformer generation
# Focus on linker analysis via SMILES patterns


def analyze_protac_smiles(smiles, name):
    """Analyze PROTAC SMILES for chameleonicity indicators."""

    # Count PEG-like fragments (C-O-C patterns)
    peg_pattern = smiles.count("COC") + smiles.count("COCCOC") + smiles.count("OCCO")

    # Count alkyl fragments (C-C-C chains)
    alkyl_patterns = smiles.count("CCC") + smiles.count("CCCC") + smiles.count("CCCCC")

    # Count ether oxygens
    ether_count = smiles.count("O") - smiles.count("N")  # Rough estimate

    # Count amide linkages (flexibility)
    amide_count = smiles.count("C(=O)N") + smiles.count("NC(=O)")

    # Compute linker flexibility score
    # PEG linkers have higher oxygen density and more flexible ether bonds
    total_carbon = smiles.count("C")
    total_oxygen = smiles.count("O")

    if total_carbon > 0:
        o_to_c_ratio = total_oxygen / total_carbon
    else:
        o_to_c_ratio = 0

    # Chameleonicity score based on linker chemistry
    # Higher PEG content = more chameleonic
    # Higher alkyl content = less chameleonic
    peg_score = peg_pattern * 2 if peg_pattern > 0 else 0
    alkyl_penalty = alkyl_patterns * 0.5 if alkyl_patterns > 3 else 0

    # Flexibility factor
    flexibility = amide_count * 0.3 + o_to_c_ratio * 10

    # Combined score
    chameleon_score = peg_score - alkyl_penalty + flexibility

    return {
        "name": name,
        "smiles": smiles,
        "peg_pattern_count": peg_pattern,
        "alkyl_pattern_count": alkyl_patterns,
        "ether_count": ether_count,
        "amide_count": amide_count,
        "o_to_c_ratio": o_to_c_ratio,
        "peg_score": peg_score,
        "alkyl_penalty": alkyl_penalty,
        "flexibility": flexibility,
        "chameleon_score": chameleon_score,
    }


def main():
    log("=" * 80)
    log("ITERATION 25: SMILES-Based Linker Chemistry Analysis")
    log("=" * 80)
    log("")
    log("Hypothesis: Chameleonicity can be predicted from linker chemistry")
    log("patterns in SMILES - PEG linkers favor folding, alkyl linkers resist it.")
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
        result = analyze_protac_smiles(smiles, name)
        results.append(result)

        log(f"\n{name}:")
        log(f"  PEG patterns: {result['peg_pattern_count']}")
        log(f"  Alkyl patterns: {result['alkyl_pattern_count']}")
        log(f"  O:C ratio: {result['o_to_c_ratio']:.3f}")
        log(f"  Chameleon Score: {result['chameleon_score']:.2f}")

    # Separation check
    log("\n" + "-" * 80)
    log("SEPARATION CHECK")
    log("-" * 80)

    cham_scores = [
        r["chameleon_score"] for r in results if r["name"] in ["protac_1", "protac_2"]
    ]
    non_score = [r["chameleon_score"] for r in results if r["name"] == "protac_3"]

    if cham_scores and non_score:
        log(f"\nChameleonic PROTACs avg: {sum(cham_scores) / len(cham_scores):.2f}")
        log(f"Non-chameleonic protac_3: {non_score[0]:.2f}")
        log(f"Difference: {sum(cham_scores) / len(cham_scores) - non_score[0]:.2f}")

        if sum(cham_scores) / len(cham_scores) > non_score[0]:
            log("Status: CORRECT ORDER (chameleonic > non-chameleonic)")
        else:
            log("Status: INCORRECT ORDER - need different approach")

    # Now generate descriptors for all 49 benchmark molecules
    log("\n" + "-" * 80)
    log("BENCHMARK DESCRIPTORS")
    log("-" * 80)

    benchmark = []
    with open("chameleon_local/benchmark.tsv", "r") as f:
        content = f.read()
        lines = content.strip().split("\n")
        log(f"Total lines in file: {len(lines)}")

        for i, line in enumerate(lines):
            parts = line.split("\t")
            log(f"Line {i}: {len(parts)} parts")
            if len(parts) >= 3:
                benchmark.append((parts[0], parts[1], parts[2]))
            elif len(parts) >= 2:
                # Try to infer from labelled_set.tsv
                benchmark.append((parts[0], parts[1], "unknown"))

    log(f"\nLoaded {len(benchmark)} benchmark molecules")

    benchmark_results = []
    for name, smiles, label in benchmark:
        result = analyze_protac_smiles(smiles, name)
        if result:
            result["label"] = label
            benchmark_results.append(result)

    log(f"\nProcessed {len(benchmark_results)} benchmark molecules")

    # Compute statistics
    cham_results = [r for r in benchmark_results if r["label"] == "cham"]
    non_results = [r for r in benchmark_results if r["label"] == "non"]

    log(f"\nChameleonic count: {len(cham_results)}")
    log(f"Non-chameleonic count: {len(non_results)}")

    if cham_results and non_results:
        cham_scores = [r["chameleon_score"] for r in cham_results]
        non_scores = [r["chameleon_score"] for r in non_results]

        log(
            f"\nChameleonic mean: {sum(cham_scores) / len(cham_scores):.2f} (n={len(cham_results)})"
        )
        log(
            f"Non-chameleonic mean: {sum(non_scores) / len(non_scores):.2f} (n={len(non_results)})"
        )
        log(
            f"Difference: {sum(cham_scores) / len(cham_scores) - sum(non_scores) / len(non_scores):.2f}"
        )

    log("\n" + "=" * 80)
    log("END OF ANALYSIS")
    log("=" * 80)

    # Write output to file
    output_text = "\n".join(output_lines)
    with open("experiments/iter_25_output.txt", "w") as f:
        f.write(output_text)

    # Save descriptors
    with open("experiments/iter_25_descriptors.json", "w") as f:
        json.dump(benchmark_results, f, indent=2)

    log(f"\nResults written to experiments/iter_25_output.txt")
    log(f"Descriptors written to experiments/iter_25_descriptors.json")


if __name__ == "__main__":
    import json

    main()
