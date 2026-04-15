import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import json
from rdkit import Chem

print("Starting dihedral preference analysis...")

# Empirical dihedral propensities from crystallographic data
DIHEDRAL_PROPENSITIES = {
    "COCC": 0.72,  # Ether - strong gauche
    "CCOC": 0.68,  # Ether - strong gauche
    "OCCO": 0.75,  # Ether-ether
    "CCNC": 0.45,  # Amide
    "CONC": 0.55,  # Amide-ether
    "NCOC": 0.50,  # Amide-ether
    "CCCC": -0.35,  # Alkyl - anti
    "CCCH": -0.30,  # Alkyl - weak anti
}


def get_atom_type(atom):
    """Get simplified atom type."""
    symbol = atom.GetSymbol()
    if symbol == "C":
        return "C"
    elif symbol == "O":
        return "O"
    elif symbol == "N":
        return "N"
    else:
        return "X"


def compute_score(smiles, name):
    """Compute dihedral preference score."""
    print(f"Processing {name}...")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  Failed to parse SMILES")
        return None

    mol = Chem.AddHs(mol)

    # Count ether and alkyl patterns in linker bonds
    ether_score = 0.0
    alkyl_score = 0.0
    ether_count = 0
    alkyl_count = 0

    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        type1 = get_atom_type(atom1)
        type2 = get_atom_type(atom2)

        # Get neighbors (excluding bond partner)
        neighbors1 = [n for n in atom1.GetNeighbors() if n.GetIdx() != atom2.GetIdx()]
        neighbors2 = [n for n in atom2.GetNeighbors() if n.GetIdx() != atom1.GetIdx()]

        # Look for dihedral patterns
        for n1 in neighbors1:
            for n2 in neighbors2:
                n1_type = get_atom_type(n1)
                n2_type = get_atom_type(n2)

                # Build 4-atom pattern
                pattern = n1_type + type1 + type2 + n2_type
                reverse = n2_type + type2 + type1 + n1_type

                # Check for ether patterns (C-O-C-C)
                if pattern in ["COCC", "CCOC", "OCCO"] or reverse in [
                    "COCC",
                    "CCOC",
                    "OCCO",
                ]:
                    ether_score += 0.70
                    ether_count += 1
                # Check for alkyl patterns (C-C-C-C)
                elif pattern == "CCCC" or reverse == "CCCC":
                    alkyl_score += -0.35
                    alkyl_count += 1

    total = ether_count + alkyl_count
    if total == 0:
        total = 1  # Avoid div by zero

    # Compute scores
    ether_ratio = ether_count / total if total > 0 else 0
    alkyl_ratio = alkyl_count / total if total > 0 else 0

    mean_score = (ether_score + alkyl_score) / total if total > 0 else 0
    composite = mean_score * 0.6 + (ether_ratio * 0.7 - alkyl_ratio * 0.3) * 0.4

    return {
        "name": name,
        "smiles": smiles,
        "mean_dihedral_pref": float(mean_score),
        "ether_dihedrals": ether_count,
        "alkyl_dihedrals": alkyl_count,
        "ether_ratio": float(ether_ratio),
        "alkyl_ratio": float(alkyl_ratio),
        "composite_score": float(composite),
    }


def load_molecules(filepath):
    """Load molecules from TSV."""
    molecules = []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split(None, 1)
            if len(parts) >= 2:
                name = parts[0]
                rest = parts[1]
                smile_end = rest.find(" ")
                if smile_end == -1:
                    smiles = rest
                else:
                    smiles = rest[:smile_end]
                molecules.append((name, smiles))
    return molecules


print("=" * 60)
print("ITERATION 26: DIHEDRAL PREFERENCE ANALYSIS (No NumPy)")
print("=" * 60)

# Load molecules
user_protacs = load_molecules("chameleon_local/user_protacs.tsv")
benchmark = load_molecules("chameleon_local/benchmark.tsv")

output_lines = []
all_descriptors = []

# Analyze user PROTACs
output_lines.append("\nUSER PROTACs ANALYSIS")
output_lines.append("-" * 60)

protac_results = []
for name, smiles in user_protacs:
    result = compute_score(smiles, name)
    if result:
        protac_results.append(result)
        output_lines.append(f"\n{name}:")
        output_lines.append(
            f"  Mean Dihedral Pref: {result['mean_dihedral_pref']:+.3f}"
        )
        output_lines.append(f"  Composite Score:    {result['composite_score']:+.3f}")
        output_lines.append(f"  Ether Dihedrals:    {result['ether_dihedrals']}")
        output_lines.append(f"  Alkyl Dihedrals:    {result['alkyl_dihedrals']}")
        output_lines.append(f"  Ether Ratio:        {result['ether_ratio']:.1%}")

# Summary
if len(protac_results) >= 3:
    output_lines.append("\n" + "-" * 60)
    output_lines.append("PROTAC COMPARISON SUMMARY")
    output_lines.append("-" * 60)

    for r in protac_results:
        label = "CHAM" if r["name"] in ["protac_1", "protac_2"] else "NON"
        output_lines.append(
            f"{r['name']:12s} | Composite: {r['composite_score']:+.3f} | Ether%: {r['ether_ratio']:.1%} | {label}"
        )

    # Check separation
    cham_scores = [
        r["composite_score"]
        for r in protac_results
        if r["name"] in ["protac_1", "protac_2"]
    ]
    non_scores = [
        r["composite_score"] for r in protac_results if r["name"] == "protac_3"
    ]

    if cham_scores and non_scores:
        cham_avg = sum(cham_scores) / len(cham_scores)
        non_avg = sum(non_scores) / len(non_scores)
        gap = cham_avg - non_avg
        output_lines.append(f"\nChameleonic avg: {cham_avg:+.3f}")
        output_lines.append(f"Non-chameleonic avg: {non_avg:+.3f}")
        output_lines.append(f"Separation gap: {gap:+.3f}")

        if gap > 0.1:
            output_lines.append("\n[SUCCESS] Clear separation achieved!")
        else:
            output_lines.append("\n[LIMITED] Gap is small")

# Analyze full benchmark
output_lines.append("\n" + "=" * 60)
output_lines.append("FULL BENCHMARK ANALYSIS")
output_lines.append("=" * 60)

for name, smiles in benchmark:
    result = compute_score(smiles, name)
    if result:
        all_descriptors.append(result)

output_lines.append(f"\nComputed descriptors for {len(all_descriptors)} molecules")

# Statistics
if all_descriptors:
    scores = [d["composite_score"] for d in all_descriptors]
    mean_score = sum(scores) / len(scores)
    min_score = min(scores)
    max_score = max(scores)
    output_lines.append(f"\nComposite Score Statistics:")
    output_lines.append(f"  Mean: {mean_score:+.3f}")
    output_lines.append(f"  Range: [{min_score:+.3f}, {max_score:+.3f}]")

# Write output
output_text = "\n".join(output_lines)

with open("experiments/iter_26_output.txt", "w") as f:
    f.write(output_text)
print(output_text)

# Save descriptors
with open("experiments/iter_26_descriptors.json", "w") as f:
    json.dump(all_descriptors, f, indent=2)

print("\n\nResults saved to experiments/iter_26_output.txt")
print("Descriptors saved to experiments/iter_26_descriptors.json")
