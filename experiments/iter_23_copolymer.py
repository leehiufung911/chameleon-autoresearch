import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"

import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

# Debug: print to stderr immediately
print("DEBUG: Script started", file=sys.stderr)


def identify_blocks(mol):
    """Identify hydrophilic and hydrophobic blocks in the molecule."""
    atoms = [atom for atom in mol.GetAtoms()]
    n_atoms = len(atoms)

    # Classify each atom
    atom_types = []
    for atom in atoms:
        symbol = atom.GetSymbol()
        is_aromatic = atom.GetIsAromatic()

        # Hydrophilic atoms: N, O (not aromatic)
        if symbol in ["N", "O"] and not is_aromatic:
            atom_types.append("hydrophilic")
        else:
            atom_types.append("hydrophobic")

    # Identify contiguous blocks
    blocks = []
    if n_atoms == 0:
        return blocks

    current_type = atom_types[0]
    current_start = 0

    for i in range(1, n_atoms):
        if atom_types[i] != current_type:
            blocks.append(
                {
                    "type": current_type,
                    "start": current_start,
                    "end": i - 1,
                    "size": i - current_start,
                }
            )
            current_type = atom_types[i]
            current_start = i

    # Add last block
    blocks.append(
        {
            "type": current_type,
            "start": current_start,
            "end": n_atoms - 1,
            "size": n_atoms - current_start,
        }
    )

    return blocks


def compute_hydrophilic_connectivity(blocks):
    """Compute how well hydrophilic blocks are connected."""
    philic_blocks = [b for b in blocks if b["type"] == "hydrophilic"]
    phobic_blocks = [b for b in blocks if b["type"] == "hydrophobic"]

    if len(philic_blocks) <= 1:
        return 1.0

    # Count hydrophilic block adjacencies
    adjacent_count = 0

    for i, block in enumerate(philic_blocks[:-1]):
        next_block = philic_blocks[i + 1]
        has_phobic_between = any(
            phobic["start"] > block["end"] and phobic["end"] < next_block["start"]
            for phobic in phobic_blocks
        )
        if not has_phobic_between:
            adjacent_count += 1

    connectivity = adjacent_count / max(1, len(philic_blocks) - 1)
    return connectivity


def compute_block_chi(blocks):
    """Compute Flory-Huggins chi-like parameter between blocks."""
    philic_blocks = [b for b in blocks if b["type"] == "hydrophilic"]
    phobic_blocks = [b for b in blocks if b["type"] == "hydrophobic"]

    if not philic_blocks or not phobic_blocks:
        return 0.0

    avg_philic_size = np.mean([b["size"] for b in philic_blocks])
    avg_phobic_size = np.mean([b["size"] for b in phobic_blocks])

    size_ratio = abs(avg_philic_size - avg_phobic_size) / max(
        avg_philic_size, avg_phobic_size
    )
    return size_ratio


def compute_copolymer_score(smiles, name):
    """Compute block copolymer chameleonicity score."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    blocks = identify_blocks(mol)

    if not blocks:
        return None

    connectivity = compute_hydrophilic_connectivity(blocks)
    chi = compute_block_chi(blocks)

    total_atoms = sum(b["size"] for b in blocks)
    philic_atoms = sum(b["size"] for b in blocks if b["type"] == "hydrophilic")
    fraction_philic = philic_atoms / total_atoms if total_atoms > 0 else 0

    n_blocks = len(blocks)

    # Chameleonicity score
    score = connectivity * (1 + chi) * fraction_philic * (1 / max(1, n_blocks / 5))

    return {
        "name": name,
        "smiles": smiles,
        "n_blocks": n_blocks,
        "philic_fraction": fraction_philic,
        "hydrophilic_connectivity": connectivity,
        "block_chi": chi,
        "copolymer_score": score,
    }


def main():
    print("Starting Block Copolymer Analysis...", file=sys.stderr)

    benchmark_file = "chameleon_local/benchmark.tsv"
    user_protacs_file = "chameleon_local/user_protacs.tsv"

    results = []

    # Process benchmark molecules
    try:
        with open(benchmark_file, "r") as f:
            lines = f.readlines()[1:]
            for line in lines:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    name = parts[0]
                    smiles = parts[1]
                    result = compute_copolymer_score(smiles, name)
                    if result:
                        results.append(result)
        print(f"Processed {len(results)} benchmark molecules", file=sys.stderr)
    except Exception as e:
        print(f"Error reading benchmark: {e}", file=sys.stderr)

    # Process user PROTACs
    try:
        with open(user_protacs_file, "r") as f:
            lines = f.readlines()[1:]
            for line in lines:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    name = parts[0]
                    smiles = parts[1]
                    result = compute_copolymer_score(smiles, name)
                    if result:
                        results.append(result)
        print(f"Total molecules: {len(results)}", file=sys.stderr)
    except Exception as e:
        print(f"Error reading user_protacs: {e}", file=sys.stderr)

    # Output results
    output_lines = []
    output_lines.append("=" * 70)
    output_lines.append("ITERATION 23: Block Copolymer Theory Analysis")
    output_lines.append("=" * 70)
    output_lines.append("")

    # Focus on PROTACs
    protacs = [r for r in results if r["name"].startswith("protac_")]
    output_lines.append(f"PROTAC Results ({len(protacs)} found):")
    output_lines.append("-" * 70)

    for protac in sorted(protacs, key=lambda x: x["name"]):
        output_lines.append(f"\n{protac['name']}:")
        output_lines.append(f"  Blocks: {protac['n_blocks']}")
        output_lines.append(f"  Hydrophilic fraction: {protac['philic_fraction']:.3f}")
        output_lines.append(
            f"  Hydrophilic connectivity: {protac['hydrophilic_connectivity']:.3f}"
        )
        output_lines.append(f"  Block chi: {protac['block_chi']:.3f}")
        output_lines.append(f"  Copolymer Score: {protac['copolymer_score']:.4f}")

    # Calculate separation
    cham_scores = [
        r["copolymer_score"] for r in protacs if r["name"] in ["protac_1", "protac_2"]
    ]
    non_scores = [r["copolymer_score"] for r in protacs if r["name"] == "protac_3"]

    if cham_scores and non_scores:
        avg_cham = np.mean(cham_scores)
        avg_non = np.mean(non_scores)
        separation = avg_cham - avg_non

        output_lines.append(f"\n{'=' * 70}")
        output_lines.append("PROTAC SEPARATION ANALYSIS")
        output_lines.append(f"{'=' * 70}")
        output_lines.append(f"Chameleonic (protac_1/2) average: {avg_cham:.4f}")
        output_lines.append(f"Non-chameleonic (protac_3) average: {avg_non:.4f}")
        output_lines.append(f"Separation: {separation:+.4f}")

        if separation > 0:
            output_lines.append("\nSUCCESS: Chameleonic PROTACs score higher")
        else:
            output_lines.append("\nFAILED: Non-chameleonic scores higher")

    output_text = "\n".join(output_lines)

    # Write to output file
    with open("experiments/iter_23_output.txt", "w") as f:
        f.write(output_text)

    print(output_text)

    # Save descriptors
    descriptors_data = []
    for r in results:
        descriptors_data.append(
            {
                "name": r["name"],
                "copolymer_score": r["copolymer_score"],
                "n_blocks": r["n_blocks"],
                "philic_fraction": r["philic_fraction"],
                "hydrophilic_connectivity": r["hydrophilic_connectivity"],
                "block_chi": r["block_chi"],
            }
        )

    with open("experiments/iter_23_descriptors.json", "w") as f:
        json.dump(descriptors_data, f, indent=2)

    print(
        "\nDescriptors saved to experiments/iter_23_descriptors.json", file=sys.stderr
    )


if __name__ == "__main__":
    main()
