import sys
import os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import json

print("=" * 70, file=sys.stderr)
print("ITERATION 28: Ramachandran-like Dihedral Preference Analysis", file=sys.stderr)
print("=" * 70, file=sys.stderr)
print("", file=sys.stderr)
print(
    "Hypothesis: Chameleonicity can be predicted from intrinsic dihedral",
    file=sys.stderr,
)
print(
    "preferences of linker bonds. PEG linkers favor gauche (folded) conformations",
    file=sys.stderr,
)
print("while alkyl linkers favor anti (extended) conformations.", file=sys.stderr)
print("", file=sys.stderr)


def extract_linker_dihedrals(mol):
    """Extract dihedral angles from conformer to create Ramachandran-like plot."""
    if mol.GetNumConformers() == 0:
        return None

    conf = mol.GetConformer(0)
    dihedrals = []

    # Find all rotatable bonds (non-ring, non-terminal)
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue

        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()

        # Get neighbors for dihedral
        atom1 = mol.GetAtomWithIdx(a1)
        atom2 = mol.GetAtomWithIdx(a2)

        n1 = [n.GetIdx() for n in atom1.GetNeighbors() if n.GetIdx() != a2]
        n2 = [n.GetIdx() for n in atom2.GetNeighbors() if n.GetIdx() != a1]

        if len(n1) > 0 and len(n2) > 0:
            try:
                # Take first neighbors to form dihedral
                dihedral = AllChem.GetDihedral(conf, n1[0], a1, a2, n2[0])

                # Classify bond type
                t1 = atom1.GetAtomicNum()
                t2 = atom2.GetAtomicNum()

                bond_type = None
                if t1 == 6 and t2 == 8:
                    bond_type = "C-O"
                elif t1 == 8 and t2 == 6:
                    bond_type = "O-C"
                elif t1 == 6 and t2 == 6:
                    bond_type = "C-C"
                elif t1 == 6 and t2 == 7:
                    bond_type = "C-N"
                elif t1 == 7 and t2 == 6:
                    bond_type = "N-C"

                if bond_type:
                    dihedrals.append(
                        {
                            "bond_type": bond_type,
                            "dihedral": float(np.degrees(dihedral)),
                        }
                    )
            except:
                pass

    return dihedrals


def compute_gauche_preference(mol):
    """Compute fraction of gauche vs anti conformations."""
    dihedrals = extract_linker_dihedrals(mol)
    if not dihedrals:
        return None

    co_dihedrals = [
        d["dihedral"] for d in dihedrals if d["bond_type"] in ["C-O", "O-C"]
    ]
    cc_dihedrals = [d["dihedral"] for d in dihedrals if d["bond_type"] == "C-C"]

    results = {"co_dihedrals": co_dihedrals, "cc_dihedrals": cc_dihedrals}

    # Classify as gauche (60±30° or 300±30°) or anti (180±30°)
    def classify_gauche(dihed):
        dihed = abs(dihed)
        if dihed < 30:  # syn
            return "syn"
        elif 30 <= dihed <= 90 or 270 <= dihed <= 330:  # gauche
            return "gauche"
        elif 150 <= dihed <= 210:  # anti
            return "anti"
        return "other"

    co_types = [classify_gauche(d) for d in co_dihedrals]
    cc_types = [classify_gauche(d) for d in cc_dihedrals]

    results["co_gauche_frac"] = (
        co_types.count("gauche") / len(co_types) if co_types else 0
    )
    results["co_anti_frac"] = co_types.count("anti") / len(co_types) if co_types else 0
    results["cc_gauche_frac"] = (
        cc_types.count("gauche") / len(cc_types) if cc_types else 0
    )
    results["cc_anti_frac"] = cc_types.count("anti") / len(cc_types) if cc_types else 0

    # Preference score: PEG favors gauche, alkyl favors anti
    # Normalize by total linker composition
    n_co = len(co_dihedrals)
    n_cc = len(cc_dihedrals)
    total = n_co + n_cc

    if total > 0:
        # Weighted preference: CO-gauche is positive, CC-anti is negative
        pref = (
            results["co_gauche_frac"] * n_co - results["cc_anti_frac"] * n_cc
        ) / total
    else:
        pref = 0

    results["preference_score"] = pref
    results["n_co"] = n_co
    results["n_cc"] = n_cc

    return results


def main():
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

    print("Processing PROTACs...", file=sys.stderr)

    results = []

    for name, smiles in protacs:
        print(f"Processing {name}...", file=sys.stderr)

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  Failed to parse SMILES", file=sys.stderr)
            continue

        mol = Chem.AddHs(mol)

        # Generate a single conformer for dihedral analysis
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 0

        AllChem.EmbedMolecule(mol, params)
        if mol.GetNumConformers() == 0:
            print(f"  Failed to embed conformer", file=sys.stderr)
            continue

        dihedral_data = compute_gauche_preference(mol)
        if dihedral_data is None:
            print(f"  No dihedrals extracted", file=sys.stderr)
            continue

        result = {
            "name": name,
            "label": labels.get(name, "?"),
            "preference_score": float(dihedral_data["preference_score"]),
            "n_co": int(dihedral_data["n_co"]),
            "n_cc": int(dihedral_data["n_cc"]),
            "co_gauche_frac": float(dihedral_data["co_gauche_frac"]),
            "cc_anti_frac": float(dihedral_data["cc_anti_frac"]),
        }
        results.append(result)

        print(
            f"  {result['label']}: Pref={result['preference_score']:.3f}, "
            f"CO={result['n_co']}, CC={result['n_cc']}",
            file=sys.stderr,
        )

    # Summary
    print("", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print("SUMMARY", file=sys.stderr)
    print("=" * 70, file=sys.stderr)

    if len(results) >= 3:
        cham_scores = [r["preference_score"] for r in results if r["label"] == "CHAM"]
        non_scores = [r["preference_score"] for r in results if r["label"] == "NON"]

        cham_mean = np.mean(cham_scores) if cham_scores else 0
        non_mean = np.mean(non_scores) if non_scores else 0
        gap = cham_mean - non_mean

        print(f"Chameleonic avg preference: {cham_mean:.3f}", file=sys.stderr)
        print(f"Non-chameleonic avg: {non_mean:.3f}", file=sys.stderr)
        print(f"Gap: {gap:+.3f}", file=sys.stderr)

        if gap > 0.05:
            print(
                "Verdict: SUCCESS - PEG-linked PROTACs favor gauche conformations",
                file=sys.stderr,
            )
        elif gap > 0:
            print("Verdict: PARTIAL SUCCESS", file=sys.stderr)
        else:
            print("Verdict: FAILED - Alkyl-linked scores higher", file=sys.stderr)

    # Write output
    output_lines = []
    output_lines.append("ITERATION 28: Ramachandran-like Dihedral Preference Analysis")
    output_lines.append("")
    output_lines.append("Results:")
    for r in results:
        output_lines.append(
            f"{r['name']} ({r['label']}): Pref={r['preference_score']:.3f}"
        )
    output_lines.append("")
    if len(results) >= 3:
        output_lines.append(f"Gap: {gap:+.3f}")

    output_text = "\n".join(output_lines)

    with open("experiments/iter_28_output.txt", "w") as f:
        f.write(output_text)

    with open("experiments/iter_28_descriptors.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\nOutput written to experiments/iter_28_output.txt", file=sys.stderr)
    print(
        "Descriptors written to experiments/iter_28_descriptors.json", file=sys.stderr
    )


if __name__ == "__main__":
    main()
