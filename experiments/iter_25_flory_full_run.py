#!/usr/bin/env python3
"""
Iteration 25: Flory-Huggins Theory Analysis for Chameleonicity

Hypothesis: Chameleonicity can be predicted from 2D topology using
Flory-Huggins-inspired interaction density. PEG linkers create more
interaction sites per flexible unit than alkyl linkers.

Physical basis:
- PEG (gauche) enables folding by creating intramolecular interaction sites
- Alkyl (anti) resists folding with fewer interaction sites
- Chi score = interaction density = total_interactions / flexibility
- Ether preference = linker_ether - linker_alkyl (positive = PEG-rich)

Novel approach: Pure 2D SMILES analysis - no 3D conformers, no implicit solvent bias
"""

import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import json

output_lines = []


def log(msg):
    output_lines.append(str(msg))
    print(str(msg), flush=True)


def analyze_smiles(smiles, name):
    """Analyze SMILES for chameleonicity indicators."""
    # Count ether oxygens (-O-)
    ether_count = smiles.count("OCC") + smiles.count("CCO") + smiles.count("COC")

    # Count amide groups (H-bond forming)
    amide_count = smiles.count("N(C") + smiles.count("C(=O)N")

    # Count polar oxygens (=O)
    polar_oxygens = smiles.count("=O")

    # Linker patterns
    linker_ether = smiles.count("COCC") + smiles.count("CCOC") + smiles.count("COCCO")
    linker_alkyl = smiles.count("CCCC") + smiles.count("CCCCC")

    # Total interaction sites
    total_interactions = ether_count + amide_count * 2 + polar_oxygens

    # Chain flexibility (rough estimate from carbon count)
    chain_flexibility = smiles.count("C") - 5

    # Flory-Huggins chi score: interaction density
    if chain_flexibility > 0:
        chi_score = total_interactions / chain_flexibility
    else:
        chi_score = total_interactions

    # Ether preference
    ether_preference = linker_ether - linker_alkyl

    # Combined chameleonicity score
    cham_score = chi_score + ether_preference * 0.5

    return {
        "name": name,
        "smiles": smiles,
        "ether_count": ether_count,
        "amide_count": amide_count,
        "polar_oxygens": polar_oxygens,
        "linker_ether": linker_ether,
        "linker_alkyl": linker_alkyl,
        "total_interactions": total_interactions,
        "flexibility": max(0, chain_flexibility),
        "chi_score": chi_score,
        "ether_preference": ether_preference,
        "cham_score": cham_score,
    }


def main():
    log("=" * 80)
    log("ITERATION 25: Flory-Huggins Theory Analysis")
    log("2D Polymer Physics Approach")
    log("=" * 80)
    log("")
    log("Hypothesis: Chameleonicity from 2D topology using")
    log("Flory-Huggins-inspired interaction density.")
    log("")

    # Load user PROTACs
    log("-" * 80)
    log("USER PROTACS ANALYSIS")
    log("-" * 80)

    protacs = {}
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                protacs[parts[0]] = parts[1]

    labels = {"protac_1": "cham", "protac_2": "cham", "protac_3": "non"}

    results = []
    for name, smiles in protacs.items():
        result = analyze_smiles(smiles, name)
        result["label"] = labels[name]
        results.append(result)

        log(f"\n{name} ({result['label']}):")
        log(f"  Ether sites: {result['ether_count']}")
        log(f"  Amide sites: {result['amide_count']}")
        log(f"  Ether preference: {result['ether_preference']}")
        log(f"  Chi score: {result['chi_score']:.3f}")
        log(f"  Cham score: {result['cham_score']:.3f}")

    # Check separation
    log("\n" + "-" * 80)
    log("SEPARATION CHECK")
    log("-" * 80)

    cham_scores = [r["cham_score"] for r in results if r["label"] == "cham"]
    non_scores = [r["cham_score"] for r in results if r["label"] == "non"]

    if cham_scores and non_scores:
        cham_mean = sum(cham_scores) / len(cham_scores)
        non_mean = sum(non_scores) / len(non_scores)
        log(f"\nChameleonic avg: {cham_mean:.3f}")
        log(f"Non-chameleonic: {non_mean:.3f}")
        log(f"Difference: {cham_mean - non_mean:.3f}")
        if cham_mean > non_mean:
            log("Status: CORRECT ORDER")
        else:
            log("Status: REVERSED ORDER")

    # Load full benchmark
    log("\n" + "-" * 80)
    log("FULL BENCHMARK ANALYSIS")
    log("-" * 80)

    benchmark = []
    with open("chameleon_local/labelled_set.tsv", "r") as f:
        lines = f.readlines()
        for line in lines[1:]:  # Skip header
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                benchmark.append((parts[0], parts[1], parts[2]))

    log(f"Loaded {len(benchmark)} molecules from labelled_set.tsv")

    # Process all benchmark molecules
    bench_results = []
    for name, smiles, label in benchmark:
        result = analyze_smiles(smiles, name)
        result["label"] = label
        bench_results.append(result)

    log(f"Processed {len(bench_results)} molecules")

    # Statistics
    import numpy as np

    cham_results = [r for r in bench_results if r["label"] == "cham"]
    non_results = [r for r in bench_results if r["label"] == "non"]

    log(f"\nChameleonic: {len(cham_results)}")
    log(f"Non-chameleonic: {len(non_results)}")

    if cham_results and non_results:
        cham_scores = [r["cham_score"] for r in cham_results]
        non_scores = [r["cham_score"] for r in non_results]

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

    # Save descriptors for all molecules
    descriptors = []
    for r in bench_results:
        descriptors.append(
            {
                "name": r["name"],
                "label": r["label"],
                "cham_score": r["cham_score"],
                "chi_score": r["chi_score"],
                "ether_preference": r["ether_preference"],
                "ether_count": r["ether_count"],
                "amide_count": r["amide_count"],
                "linker_ether": r["linker_ether"],
                "linker_alkyl": r["linker_alkyl"],
            }
        )

    with open("experiments/iter_25_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)

    log(f"\nResults: experiments/iter_25_output.txt")
    log(f"Descriptors: experiments/iter_25_descriptors.json")


if __name__ == "__main__":
    main()
