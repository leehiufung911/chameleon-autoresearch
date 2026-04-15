#!/usr/bin/env python
import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys

sys.path.insert(
    0, "C:/Users/mic23/prototype-embed/chameleon-research-loop/chameleon_local"
)

from rdkit import Chem
import numpy as np
import json

print("=" * 80, flush=True)
print("ITERATION 22: Fragment-Based Chameleonicity Prediction", flush=True)
print("=" * 80, flush=True)
print("", flush=True)

# SMARTS patterns
PROMOTING = {
    "peg_chain": "[OX2][C;!R][C;!R][OX2]",
    "ether": "[OX2]-[!H]",
    "amide": "[NX3](C=O)",
    "ester": "C(=O)[OX2]",
}

INHIBITING = {
    "long_alkyl": "[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]",
    "longer_alkyl": "[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]",
}


def count_matches(mol, patterns):
    counts = {}
    for name, smarts in patterns.items():
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat:
                counts[name] = len(mol.GetSubstructMatches(pat))
            else:
                counts[name] = 0
        except Exception as e:
            counts[name] = 0
    return counts


def score_mol(mol):
    prom = count_matches(mol, PROMOTING)
    inhib = count_matches(mol, INHIBITING)

    prom_score = (
        prom.get("peg_chain", 0) * 3.0
        + prom.get("amide", 0) * 2.0
        + prom.get("ether", 0) * 1.0
        + prom.get("ester", 0) * 0.5
    )

    inhib_score = inhib.get("long_alkyl", 0) * 2.0 + inhib.get("longer_alkyl", 0) * 3.0

    total = prom_score + inhib_score
    balance = prom_score / total if total > 0 else 0.5

    score = balance * (1 + np.log1p(prom_score)) / (1 + np.log1p(inhib_score + 1))

    return {
        "promoting_score": prom_score,
        "inhibiting_score": inhib_score,
        "balance": balance,
        "chameleon_score": float(score),
        "prom_details": prom,
        "inhib_details": inhib,
    }


# Load PROTACs
print("Loading user PROTACs...", flush=True)
protacs = []
with open(
    "C:/Users/mic23/prototype-embed/chameleon-research-loop/chameleon_local/user_protacs.tsv",
    "r",
) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            name = parts[0]
            smiles = parts[1]
            label = "CHAM" if name in ["protac_1", "protac_2"] else "NON"
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                protacs.append({"name": name, "mol": mol, "label": label})

print(f"Loaded {len(protacs)} PROTACs", flush=True)

# Load benchmark
print("Loading benchmark molecules...", flush=True)
bench_mols = []
with open(
    "C:/Users/mic23/prototype-embed/chameleon-research-loop/chameleon_local/benchmark.tsv",
    "r",
) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            name, smiles, label = parts[0], parts[1], parts[2]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                bench_mols.append({"name": name, "mol": mol, "label": label})

print(f"Loaded {len(bench_mols)} benchmark molecules", flush=True)
print("", flush=True)

# Process all
output_lines = []
all_results = []

output_lines.append("USER PROTAC RESULTS:")
output_lines.append("-" * 80)
output_lines.append(
    f"{'Name':<15} {'Label':<8} {'Promo':>8} {'Inhib':>8} {'Balance':>8} {'Score':>10}"
)
output_lines.append("-" * 80)

for item in protacs:
    result = score_mol(item["mol"])
    result["name"] = item["name"]
    result["label"] = item["label"]
    all_results.append(result)

    output_lines.append(
        f"{item['name']:<15} {item['label']:<8} {result['promoting_score']:>8.2f} "
        f"{result['inhibiting_score']:>8.2f} {result['balance']:>8.3f} {result['chameleon_score']:>10.4f}"
    )
    output_lines.append(
        f"  PEG={result['prom_details']['peg_chain']}, "
        f"Amide={result['prom_details']['amide']}, "
        f"Ether={result['prom_details']['ether']}, "
        f"Alkyl5={result['inhib_details']['long_alkyl']}"
    )

# Process benchmark
for item in bench_mols:
    result = score_mol(item["mol"])
    result["name"] = item["name"]
    result["label"] = item["label"]
    all_results.append(result)

# Summary
cham_scores = [r["chameleon_score"] for r in all_results if r["label"] == "CHAM"]
non_scores = [r["chameleon_score"] for r in all_results if r["label"] == "NON"]

output_lines.append("")
output_lines.append("=" * 80)
output_lines.append("SUMMARY STATISTICS")
output_lines.append("=" * 80)
output_lines.append("")

if cham_scores:
    output_lines.append(
        f"Chameleonic (n={len(cham_scores)}): mean={np.mean(cham_scores):.4f}, std={np.std(cham_scores):.4f}"
    )

if non_scores:
    output_lines.append(
        f"Non-chameleonic (n={len(non_scores)}): mean={np.mean(non_scores):.4f}, std={np.std(non_scores):.4f}"
    )

# PROTAC check
protac_results = {r["name"]: r for r in all_results if r["name"].startswith("protac_")}
if len(protac_results) >= 3:
    output_lines.append("")
    output_lines.append("PROTAC SEPARATION:")
    for name in ["protac_1", "protac_2", "protac_3"]:
        if name in protac_results:
            r = protac_results[name]
            output_lines.append(
                f"  {name}: score={r['chameleon_score']:.4f} ({r['label']})"
            )

    if all(n in protac_results for n in ["protac_1", "protac_2", "protac_3"]):
        cham_avg = np.mean(
            [
                protac_results["protac_1"]["chameleon_score"],
                protac_results["protac_2"]["chameleon_score"],
            ]
        )
        non_score = protac_results["protac_3"]["chameleon_score"]
        output_lines.append(f"\n  Chameleonic avg: {cham_avg:.4f}")
        output_lines.append(f"  Non-chameleonic: {non_score:.4f}")
        output_lines.append(f"  Difference:      {cham_avg - non_score:+.4f}")
        if cham_avg > non_score:
            output_lines.append("  VERDICT: SUCCESS")
        else:
            output_lines.append("  VERDICT: FAIL")

output_lines.append("")
output_lines.append("=" * 80)
output_lines.append("Experiment complete.")

# Print and save
output_text = "\n".join(output_lines)
print(output_text, flush=True)

# Save to file
filepath = "C:/Users/mic23/prototype-embed/chameleon-research-loop/experiments/iter_22_output.txt"
with open(filepath, "w") as f:
    f.write(output_text)
print(f"\nSaved to {filepath}", flush=True)

# Save descriptors
filepath_json = "C:/Users/mic23/prototype-embed/chameleon-research-loop/experiments/iter_22_descriptors.json"
descriptors_clean = []
for r in all_results:
    clean = {
        k: float(v) if isinstance(v, (np.floating, np.integer)) else v
        for k, v in r.items()
        if k not in ["mol"]
    }
    descriptors_clean.append(clean)

with open(filepath_json, "w") as f:
    json.dump(descriptors_clean, f, indent=2)
print(f"Descriptors saved to {filepath_json}", flush=True)
