"""Iteration 27: Solvent Accessible Surface Area Entropy Analysis

Hypothesis: Chameleonic molecules show bimodal SASA distributions (high entropy)
between polar and apolar ensembles, while non-chameleonic flexible molecules
show unimodal distributions (low entropy). This captures the entropic cost
of conformational switching.
"""

import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdFreeSASA
import json
import warnings

warnings.filterwarnings("ignore")


def compute_sasa_entropy(smiles, name, n_confs=50):
    """Generate conformers and compute SASA distribution statistics

    Returns entropy metrics that capture bimodality of SASA distribution
    """
    print(f"\n{'=' * 60}")
    print(f"Processing {name}")
    print(f"{'=' * 60}")

    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Generate conformers
        print(f"Generating {n_confs} conformers...")
        AllChem.EmbedMultipleConfs(
            mol, numConfs=n_confs, randomSeed=42, params=AllChem.ETKDGv3()
        )

        if mol.GetNumConformers() == 0:
            print("WARNING: No conformers generated")
            return None

        # Compute SASA for each conformer
        print(f"Computing SASA for {mol.GetNumConformers()} conformers...")
        sasas = []
        for conf_id in range(mol.GetNumConformers()):
            sasa = rdFreeSASA.CalcSASA(mol, confIdx=conf_id)
            sasas.append(sasa)

        sasas = np.array(sasas)

        # Statistical measures
        sasa_mean = np.mean(sasas)
        sasa_std = np.std(sasas)
        sasa_range = np.max(sasas) - np.min(sasas)
        sasa_cv = sasa_std / sasa_mean if sasa_mean > 0 else 0

        # Distribution entropy (Shannon entropy of binned distribution)
        # Bin SASA values and compute entropy
        hist, bin_edges = np.histogram(sasas, bins=10, density=True)
        hist = hist[hist > 0]  # Remove zero bins
        entropy = -np.sum(hist * np.log(hist)) if len(hist) > 0 else 0

        # Bimodality coefficient
        # BC = (skew^2 + 1) / (kurtosis + 3)
        # BC > 0.55 suggests bimodality
        from scipy import stats

        skew = stats.skew(sasas)
        kurt = stats.kurtosis(sasas, fisher=False)
        bimodality = (skew**2 + 1) / (kurt + 3) if kurt > -3 else 0

        # Separation metric: ratio of 90th to 10th percentile
        separation = (
            np.percentile(sasas, 90) / np.percentile(sasas, 10)
            if np.percentile(sasas, 10) > 0
            else 1
        )

        result = {
            "name": name,
            "n_confs": len(sasas),
            "sasa_mean": float(sasa_mean),
            "sasa_std": float(sasa_std),
            "sasa_range": float(sasa_range),
            "sasa_cv": float(sasa_cv),
            "sasa_entropy": float(entropy),
            "bimodality": float(bimodality),
            "separation_ratio": float(separation),
            "sasa_min": float(np.min(sasas)),
            "sasa_max": float(np.max(sasas)),
            "sasa_values": [float(x) for x in sasas],
        }

        print(f"  SASA mean: {sasa_mean:.1f} Å²")
        print(f"  SASA std:  {sasa_std:.1f} Å²")
        print(f"  SASA CV:   {sasa_cv:.3f}")
        print(f"  Entropy:   {entropy:.3f}")
        print(f"  Bimodality: {bimodality:.3f}")

        return result

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback

        traceback.print_exc()
        return None


def main():
    print("=" * 60)
    print("Iteration 27: SASA Entropy Analysis")
    print("=" * 60)
    print()
    print("Hypothesis: Chameleonic molecules have bimodal SASA distributions")
    print("(high entropy/coefficient of variation), while non-chameleonic")
    print("flexible molecules have unimodal distributions (low entropy).")
    print()

    # Load user PROTACs
    df = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t")
    print(f"Loaded {len(df)} PROTACs")
    print()

    results = []
    descriptors_49 = []

    # Process user PROTACs
    for _, row in df.iterrows():
        result = compute_sasa_entropy(row["smiles"], row["name"], n_confs=50)
        if result:
            result["label"] = row["chameleonic"]
            results.append(result)

    # Summary table
    print("\n" + "=" * 60)
    print("SUMMARY - USER PROTACS")
    print("=" * 60)
    print()
    print(
        f"{'Name':<12} {'Label':<12} {'SASA Mean':<12} {'CV':<10} {'Entropy':<10} {'Bimodality'}"
    )
    print("-" * 70)

    for r in results:
        label_str = "CHAM" if r["label"] == "yes" else "NON"
        print(
            f"{r['name']:<12} {label_str:<12} {r['sasa_mean']:<12.1f} {r['sasa_cv']:<10.3f} {r['sasa_entropy']:<10.3f} {r['bimodality']:.3f}"
        )

    # Analysis
    print("\n" + "=" * 60)
    print("ANALYSIS")
    print("=" * 60)
    print()

    cham = [r for r in results if r["label"] == "yes"]
    non = [r for r in results if r["label"] == "no"]

    if cham and non:
        print("Chameleonic PROTACs (protac_1, protac_2):")
        for r in cham:
            print(
                f"  {r['name']}: CV={r['sasa_cv']:.3f}, Entropy={r['sasa_entropy']:.3f}, Bimodality={r['bimodality']:.3f}"
            )

        print()
        print("Non-chameleonic PROTAC (protac_3):")
        for r in non:
            print(
                f"  {r['name']}: CV={r['sasa_cv']:.3f}, Entropy={r['sasa_entropy']:.3f}, Bimodality={r['bimodality']:.3f}"
            )

        print()

        # Compare key metrics
        cham_cv = np.mean([r["sasa_cv"] for r in cham])
        non_cv = np.mean([r["sasa_cv"] for r in non])

        cham_ent = np.mean([r["sasa_entropy"] for r in cham])
        non_ent = np.mean([r["sasa_entropy"] for r in non])

        cham_bimod = np.mean([r["bimodality"] for r in cham])
        non_bimod = np.mean([r["bimodality"] for r in non])

        print(f"Average Coefficient of Variation:")
        print(f"  Chameleonic: {cham_cv:.3f}")
        print(f"  Non-chameleonic: {non_cv:.3f}")
        print(
            f"  Difference: {cham_cv - non_cv:+.3f} ({'PASS' if cham_cv > non_cv else 'FAIL'})"
        )
        print()

        print(f"Average Entropy:")
        print(f"  Chameleonic: {cham_ent:.3f}")
        print(f"  Non-chameleonic: {non_ent:.3f}")
        print(
            f"  Difference: {cham_ent - non_ent:+.3f} ({'PASS' if cham_ent > non_ent else 'FAIL'})"
        )
        print()

        print(f"Average Bimodality:")
        print(f"  Chameleonic: {cham_bimod:.3f}")
        print(f"  Non-chameleonic: {non_bimod:.3f}")
        print(f"  Difference: {cham_bimod - non_bimod:+.3f}")
        print()

        # Check for bimodality
        print("Bimodality interpretation:")
        print("  BC > 0.55 suggests bimodal distribution (two-state switching)")
        print("  BC < 0.55 suggests unimodal distribution (single state)")
        print()

        for r in results:
            bimod_status = "BIMODAL" if r["bimodality"] > 0.55 else "UNIMODAL"
            print(f"  {r['name']}: BC={r['bimodality']:.3f} -> {bimod_status}")

    # Generate descriptors for all 49 molecules
    print("\n" + "=" * 60)
    print("Generating descriptors for 49-molecule benchmark...")
    print("=" * 60)
    print()

    df_49 = pd.read_csv("chameleon_local/benchmark.tsv", sep="\t")
    print(f"Loaded {len(df_49)} molecules from benchmark")
    print()

    for idx, row in df_49.iterrows():
        if idx % 10 == 0:
            print(f"  Processing {idx + 1}/{len(df_49)}...")

        result = compute_sasa_entropy(row["smiles"], row["name"], n_confs=30)
        if result:
            descriptor = {
                "name": row["name"],
                "sasa_mean": result["sasa_mean"],
                "sasa_cv": result["sasa_cv"],
                "sasa_entropy": result["sasa_entropy"],
                "bimodality": result["bimodality"],
                "separation_ratio": result["separation_ratio"],
            }
            descriptors_49.append(descriptor)

    print(f"\nGenerated descriptors for {len(descriptors_49)} molecules")

    # Save results
    output_data = {
        "method": "SASA_Entropy_Analysis",
        "description": "SASA distribution entropy and bimodality as chameleonicity descriptors",
        "user_protacs": results,
        "benchmark_descriptors": descriptors_49,
    }

    with open("experiments/iter_27_descriptors.json", "w") as f:
        json.dump(output_data, f, indent=2)

    print("\nSaved descriptors to experiments/iter_27_descriptors.json")

    # Write output text
    output_text = f"""Iteration 27: SASA Entropy Analysis Results
{"=" * 60}

Hypothesis: Chameleonic molecules have bimodal SASA distributions (high entropy/CV),
while non-chameleonic flexible molecules have unimodal distributions (low entropy).

User PROTAC Results:
{"Name":<12} {"Label":<12} {"SASA Mean":<12} {"CV":<10} {"Entropy":<10} {"Bimodality":<10}
{"-" * 70}
"""

    for r in results:
        label_str = "CHAM" if r["label"] == "yes" else "NON"
        output_text += f"{r['name']:<12} {label_str:<12} {r['sasa_mean']:<12.1f} {r['sasa_cv']:<10.3f} {r['sasa_entropy']:<10.3f} {r['bimodality']:<10.3f}\n"

    output_text += f"""

Analysis:
Coefficient of Variation (CV = std/mean):
  Chameleonic avg: {cham_cv:.3f}
  Non-chameleonic avg: {non_cv:.3f}
  Difference: {cham_cv - non_cv:+.3f} ({"PASS - Chameleonic has higher CV" if cham_cv > non_cv else "FAIL"})

Shannon Entropy (of binned SASA distribution):
  Chameleonic avg: {cham_ent:.3f}
  Non-chameleonic avg: {non_ent:.3f}
  Difference: {cham_ent - non_ent:+.3f} ({"PASS - Chameleonic has higher entropy" if cham_ent > non_ent else "FAIL"})

Bimodality Coefficient (BC > 0.55 = bimodal):
  Chameleonic avg: {cham_bimod:.3f}
  Non-chameleonic avg: {non_bimod:.3f}
  Difference: {cham_bimod - non_bimod:+.3f}

Interpretation:
- Bimodality > 0.55 suggests two-state switching (compact ↔ extended)
- Bimodality < 0.55 suggests single dominant state
- Chameleonic molecules should show higher bimodality and entropy

Conclusion:
{"PASS: Chameleonic PROTACs show higher SASA entropy/CV than non-chameleonic" if cham_cv > non_cv and cham_ent > non_ent else "MIXED: Some metrics separate, others do not" if (cham_cv > non_cv or cham_ent > non_ent) else "FAIL: SASA entropy does not separate chameleonic from non-chameleonic PROTACs"}

Generated descriptors for {len(descriptors_49)} benchmark molecules.
Saved to experiments/iter_27_descriptors.json
"""

    with open("experiments/iter_27_output.txt", "w") as f:
        f.write(output_text)

    print(output_text)


if __name__ == "__main__":
    main()
