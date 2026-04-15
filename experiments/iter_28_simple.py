import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

from rdkit import Chem
from rdkit.Chem import Descriptors
import json


def count_hbond_sites(smiles, name):
    """Count H-bond donors and acceptors from SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    donors = 0
    acceptors = 0

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()

        # H-bond donors: N or O with H attached
        if atomic_num in [7, 8]:
            has_h = any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors())
            if has_h:
                donors += 1
            # Acceptors: N or O with lone pairs
            if len(atom.GetNeighbors()) < 4:
                acceptors += 1

    # Compute molecular properties
    mw = Descriptors.MolWt(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    return {
        "name": name,
        "mw": mw,
        "rotatable_bonds": rotatable_bonds,
        "hbd_rdkit": hbd,
        "hba_rdkit": hba,
        "hbd_manual": donors,
        "hba_manual": acceptors,
        "hb_sites_total": donors + acceptors,
        "hb_density": (donors + acceptors) / mw * 100 if mw > 0 else 0,
    }


def main():
    output_lines = []
    output_lines.append("=" * 70)
    output_lines.append("ITERATION 28: H-Bond Site Analysis (Simplified)")
    output_lines.append("=" * 70)
    output_lines.append("")

    # Load PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) >= 2:
                name = parts[0]
                smiles = parts[1].strip()
                protacs.append((name, smiles))

    labels = {"protac_1": "CHAM", "protac_2": "CHAM", "protac_3": "NON"}

    results = []
    output_lines.append("Processing PROTACs...")
    output_lines.append("-" * 70)

    for name, smiles in protacs:
        result = count_hbond_sites(smiles, name)
        if result:
            results.append(result)
            label = labels.get(name, "?")
            output_lines.append(f"\n{name} ({label}):")
            output_lines.append(f"  MW: {result['mw']:.1f}")
            output_lines.append(f"  Rotatable bonds: {result['rotatable_bonds']}")
            output_lines.append(
                f"  HBD (RDKit): {result['hbd_rdkit']}, HBA: {result['hba_rdkit']}"
            )
            output_lines.append(
                f"  HBD (manual): {result['hbd_manual']}, HBA: {result['hba_manual']}"
            )
            output_lines.append(f"  Total HB sites: {result['hb_sites_total']}")
            output_lines.append(f"  HB density: {result['hb_density']:.4f}")

    # Summary
    output_lines.append("\n" + "=" * 70)
    output_lines.append("SUMMARY")
    output_lines.append("=" * 70)

    if len(results) >= 3:
        cham_density = [
            r["hb_density"] for r in results if labels.get(r["name"]) == "CHAM"
        ]
        non_density = [
            r["hb_density"] for r in results if labels.get(r["name"]) == "NON"
        ]

        if cham_density and non_density:
            cham_mean = sum(cham_density) / len(cham_density)
            non_mean = sum(non_density) / len(non_density)
            gap = cham_mean - non_mean

            output_lines.append(f"\nChameleonic avg HB density: {cham_mean:.4f}")
            output_lines.append(f"Non-chameleonic avg: {non_mean:.4f}")
            output_lines.append(f"Separation: {gap:+.4f}")

    output_text = "\n".join(output_lines)

    # Write output
    with open("experiments/iter_28_output.txt", "w") as f:
        f.write(output_text)

    # Write descriptors
    with open("experiments/iter_28_descriptors.json", "w") as f:
        json.dump(results, f, indent=2)

    print(output_text)
    print("\nOutput written to experiments/iter_28_output.txt")
    print("Descriptors written to experiments/iter_28_descriptors.json")


if __name__ == "__main__":
    main()
