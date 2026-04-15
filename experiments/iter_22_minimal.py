import os

os.environ["PYTHONUNBUFFERED"] = "1"

output = []
output.append("=" * 80)
output.append("ITERATION 22: Fragment-Based Chameleonicity Prediction")
output.append("=" * 80)
output.append("")

import sys

sys.path.insert(
    0, "C:/Users/mic23/prototype-embed/chameleon-research-loop/chameleon_local"
)

from rdkit import Chem

output.append("Loading PROTACs...")

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

output.append(f"Loaded {len(protacs)} PROTACs")

# SMARTS patterns
PROMOTING = {
    "peg_chain": "[OX2][C;!R][C;!R][OX2]",
    "ether": "[OX2]-[!H]",
    "amide": "[NX3](C=O)",
}

INHIBITING = {
    "long_alkyl": "[C;!R]-[C;!R]-[C;!R]-[C;!R]-[C;!R]",
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
        except:
            counts[name] = 0
    return counts


output.append("")
output.append("PROTAC FRAGMENT ANALYSIS:")
output.append("-" * 80)

all_results = []

for item in protacs:
    prom = count_matches(item["mol"], PROMOTING)
    inhib = count_matches(item["mol"], INHIBITING)

    prom_score = (
        prom.get("peg_chain", 0) * 3.0
        + prom.get("amide", 0) * 2.0
        + prom.get("ether", 0) * 1.0
    )
    inhib_score = inhib.get("long_alkyl", 0) * 2.0

    total = prom_score + inhib_score
    balance = prom_score / total if total > 0 else 0.5

    import numpy as np

    score = balance * (1 + np.log1p(prom_score)) / (1 + np.log1p(inhib_score + 1))

    result = {
        "name": item["name"],
        "label": item["label"],
        "promoting_score": prom_score,
        "inhibiting_score": inhib_score,
        "balance": balance,
        "chameleon_score": float(score),
        "prom_details": prom,
        "inhib_details": inhib,
    }
    all_results.append(result)

    output.append(
        f"{item['name']}: Promo={prom_score:.1f}, Inhib={inhib_score:.1f}, Balance={balance:.3f}, Score={score:.4f}"
    )
    output.append(
        f"  PEG={prom['peg_chain']}, Amide={prom['amide']}, Ether={prom['ether']}, Alkyl5={inhib['long_alkyl']}"
    )

# Check separation
if len(all_results) >= 3:
    protac_results = {r["name"]: r for r in all_results}
    if (
        "protac_1" in protac_results
        and "protac_2" in protac_results
        and "protac_3" in protac_results
    ):
        cham_avg = (
            protac_results["protac_1"]["chameleon_score"]
            + protac_results["protac_2"]["chameleon_score"]
        ) / 2
        non_score = protac_results["protac_3"]["chameleon_score"]
        output.append("")
        output.append("PROTAC SEPARATION:")
        output.append(f"  Chameleonic avg (1+2): {cham_avg:.4f}")
        output.append(f"  Non-chameleonic (3):   {non_score:.4f}")
        output.append(f"  Difference:            {cham_avg - non_score:+.4f}")
        if cham_avg > non_score:
            output.append("  VERDICT: SUCCESS")
        else:
            output.append("  VERDICT: FAIL")

output.append("")
output.append("=" * 80)
output.append("Experiment complete.")

# Write output
output_text = "\n".join(output)
with open(
    "C:/Users/mic23/prototype-embed/chameleon-research-loop/experiments/iter_22_output.txt",
    "w",
) as f:
    f.write(output_text)

# Also print
print(output_text, flush=True)
