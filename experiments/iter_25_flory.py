import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"

# Use only basic modules
log = open("experiments/iter_25_flory.log", "w")
log.write("Step 1: Start Flory-Huggins Analysis\n")
log.flush()

log.write("Step 2: Loading PROTAC data...\n")
log.flush()

# PROTAC SMILES from the TSV file
protac_data = {
    "protac_1": "O=C(C(N1C(C2=CC=CC(NCCOCCOCCOCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=C2C1=O)=O)CC7)NC7=O",
    "protac_2": "O=C(C(N1C(C2=CC=CC(NC(COCCOCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
    "protac_3": "O=C(C(N1C(C2=CC=CC(OCC(NCCCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
}

labels = {"protac_1": "CHAM", "protac_2": "CHAM", "protac_3": "NON"}

log.write(f"Step 3: Loaded {len(protac_data)} PROTACs\n")
log.flush()

log.write("Step 4: Computing 2D descriptors...\n")
log.flush()


# Simple pattern matching on SMILES for linker analysis
def analyze_smiles(smiles):
    """Analyze SMILES for chameleonicity indicators."""
    # Count ether oxygens (-O-) vs alkyl chains
    ether_count = smiles.count("OCC") + smiles.count("CCO")

    # Count amide groups (H-bond forming)
    amide_count = smiles.count("N(C") + smiles.count("C(=O)N")

    # Count ester/amide oxygens
    polar_oxygens = smiles.count("=O")

    # Simple linker length approximation
    linker_ether = smiles.count("COCC") + smiles.count("CCOC")
    linker_alkyl = smiles.count("CCCC") + smiles.count("CCC(C)")

    # Flory-Huggins inspired: ratio of interaction sites
    # More interaction sites = more tendency to fold
    total_interactions = ether_count + amide_count * 2 + polar_oxygens

    # Estimated chain flexibility
    chain_flexibility = smiles.count("C") - 5  # rough estimate

    return {
        "ether_sites": ether_count,
        "amide_sites": amide_count,
        "polar_oxygens": polar_oxygens,
        "linker_ether": linker_ether,
        "linker_alkyl": linker_alkyl,
        "total_interactions": total_interactions,
        "flexibility": max(0, chain_flexibility),
    }


# Analyze each PROTAC
results = []
for name, smiles in protac_data.items():
    log.write(f"\nAnalyzing {name}:\n")
    log.flush()

    descriptors = analyze_smiles(smiles)
    descriptors["name"] = name
    descriptors["label"] = labels[name]

    for key, val in descriptors.items():
        if key not in ["name", "label"]:
            log.write(f"  {key}: {val}\n")
    log.flush()

    results.append(descriptors)

log.write("\nStep 5: Computing Flory-Huggins-like score...\n")
log.flush()

# Flory-Huggins inspired score
# Chi parameter analog: interaction energy per flexible unit
for r in results:
    # More interactions = more folding propensity
    # Normalized by flexibility to avoid bias toward large molecules
    if r["flexibility"] > 0:
        interaction_density = r["total_interactions"] / r["flexibility"]
        r["chi_score"] = interaction_density
    else:
        r["chi_score"] = r["total_interactions"]

    # Ether preference: positive for PEG, negative for alkyl
    r["ether_preference"] = r["linker_ether"] - r["linker_alkyl"]

    # Combined chameleonicity score
    r["cham_score"] = r["chi_score"] + r["ether_preference"] * 0.5

log.write("\nStep 6: Results Summary\n")
log.write("=" * 60 + "\n")
log.flush()

log.write(
    f"{'Molecule':<12} {'Label':<6} {'Chi':<8} {'EtherPref':<10} {'ChamScore':<10}\n"
)
log.write("-" * 60 + "\n")
log.flush()

for r in results:
    log.write(
        f"{r['name']:<12} {r['label']:<6} {r['chi_score']:<8.2f} {r['ether_preference']:<10} {r['cham_score']:<10.2f}\n"
    )
    log.flush()

log.write("\n")
log.write("Physical Interpretation:\n")
log.write("  Chi Score: Interaction density (higher = more folding sites)\n")
log.write("  Ether Preference: PEG vs alkyl linker (positive = PEG)\n")
log.write("  Cham Score: Combined (higher = more chameleonic)\n")
log.write("\n")
log.flush()

# Check separation
cham_scores = [r["cham_score"] for r in results if r["label"] == "CHAM"]
non_scores = [r["cham_score"] for r in results if r["label"] == "NON"]

if cham_scores and non_scores:
    log.write(f"Chameleonic mean: {sum(cham_scores) / len(cham_scores):.2f}\n")
    log.write(f"Non-chameleonic mean: {sum(non_scores) / len(non_scores):.2f}\n")
    log.write(
        f"Separation: {sum(cham_scores) / len(cham_scores) - sum(non_scores) / len(non_scores):.2f}\n"
    )
log.flush()

log.write("\nStep 7: Writing output files...\n")
log.flush()

# Write output
output_text = []
output_text.append("=" * 60)
output_text.append("ITERATION 25: Flory-Huggins Theory Analysis")
output_text.append("2D Polymer Physics Approach")
output_text.append("=" * 60)
output_text.append("")
output_text.append("Hypothesis: Chameleonicity can be predicted from 2D topology")
output_text.append("using Flory-Huggins-inspired interaction density.")
output_text.append("")

for r in results:
    output_text.append(f"\n{r['name']} ({r['label']}):")
    output_text.append(f"  Ether sites: {r['ether_sites']}")
    output_text.append(f"  Amide sites: {r['amide_sites']}")
    output_text.append(f"  Ether preference: {r['ether_preference']}")
    output_text.append(f"  Chi score: {r['chi_score']:.2f}")
    output_text.append(f"  Cham score: {r['cham_score']:.2f}")

output_text.append("\n" + "=" * 60)
output_text.append("Summary:")
output_text.append("=" * 60)
for r in results:
    output_text.append(f"{r['name']}: {r['cham_score']:.2f}")

with open("experiments/iter_25_output.txt", "w") as f:
    f.write("\n".join(output_text))

log.write("Output written to iter_25_output.txt\n")
log.flush()

# Write descriptors JSON
import json

with open("experiments/iter_25_descriptors.json", "w") as f:
    json.dump(results, f, indent=2)

log.write("Descriptors written to iter_25_descriptors.json\n")
log.flush()

log.write("\nDone!\n")
log.close()

print("Experiment completed", file=sys.stderr)
