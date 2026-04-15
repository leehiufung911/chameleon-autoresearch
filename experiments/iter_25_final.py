import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

from rdkit import Chem
from rdkit.Chem import AllChem
import json


def compute_rg(positions):
    import numpy as np

    positions = np.array(positions)
    com = np.mean(positions, axis=0)
    r_sq = np.sum((positions - com) ** 2, axis=1)
    rg_sq = np.mean(r_sq)
    return np.sqrt(rg_sq)


def get_conformer_energies_and_rg(mol, n_confs=15):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, n_confs, AllChem.ETKDGv3())
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")

    results = []
    for i in range(mol.GetNumConformers()):
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=i)
        ff.Initialize()
        ff.Minimize(maxIts=200)
        energy = ff.CalcEnergy()

        conf = mol.GetConformer(i)
        import numpy as np

        positions = np.array(
            [conf.GetAtomPosition(j) for j in range(mol.GetNumAtoms())]
        )
        rg = compute_rg(positions)

        results.append({"energy": energy, "rg": rg, "conf_id": i})
    return results


def analyze_protac(smiles, name):
    lines = [f"\n{'=' * 60}", f"Analyzing: {name}", f"{'=' * 60}"]

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"name": name, "status": "failed"}, "\n".join(lines)

        results = get_conformer_energies_and_rg(mol, n_confs=15)

        import numpy as np

        energies = np.array([r["energy"] for r in results])
        rg_values = np.array([r["rg"] for r in results])

        lines.append(f"\nN={len(results)} conformers")
        lines.append(
            f"Energy range: {np.min(energies):.1f} - {np.max(energies):.1f} kcal/mol"
        )
        lines.append(f"Rg range: {np.min(rg_values):.2f} - {np.max(rg_values):.2f} A")
        lines.append(f"Mean Rg: {np.mean(rg_values):.2f} A")

        if len(energies) > 2:
            corr = np.corrcoef(energies, rg_values)[0, 1]
            lines.append(f"Energy-Rg correlation: {corr:.3f}")

        return {
            "name": name,
            "status": "success",
            "mean_rg": float(np.mean(rg_values)),
            "min_rg": float(np.min(rg_values)),
            "max_rg": float(np.max(rg_values)),
            "energy_rg_corr": float(corr) if len(energies) > 2 else 0.0,
        }, "\n".join(lines)

    except Exception as e:
        return {"name": name, "status": "failed", "error": str(e)}, f"\n".join(
            lines
        ) + f"\nError: {e}"


def main():
    output = []
    output.append("=" * 60)
    output.append("ITER 25: MMFF Energy-Rg Analysis")
    output.append("=" * 60)

    # Load PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                protacs.append((parts[0], parts[1]))

    results = []
    for name, smiles in protacs:
        result, log = analyze_protac(smiles, name)
        output.append(log)
        results.append(result)

    # Summary
    output.append("\n" + "=" * 60)
    output.append("SUMMARY")
    output.append("=" * 60)
    output.append(
        f"{'Molecule':<12} {'Mean Rg':<10} {'Rg Range':<18} {'E-Rg Corr':<12}"
    )
    output.append("-" * 60)

    labels = {"protac_1": "CHAM", "protac_2": "CHAM", "protac_3": "NON"}

    for r in results:
        if r.get("status") == "success":
            name = r["name"]
            label = labels.get(name, "UNK")
            mean_rg = r["mean_rg"]
            rg_range = f"{r['min_rg']:.1f}-{r['max_rg']:.1f}"
            corr_str = f"{r.get('energy_rg_corr', 0):.3f}"
            output.append(f"{name:<12} {mean_rg:<10.2f} {rg_range:<18} {corr_str:<12}")

    output.append("\n" + "=" * 60)
    output.append("Interpretation:")
    output.append("Negative E-Rg correlation: MMFF favors compact")
    output.append("Near-zero correlation: MMFF unbiased")
    output.append("=" * 60)

    final_text = "\n".join(output)

    # Write files
    with open("experiments/iter_25_output.txt", "w") as f:
        f.write(final_text)

    descriptors = []
    for r in results:
        if r.get("status") == "success":
            descriptors.append(
                {
                    "name": r["name"],
                    "mean_rg_mmff": r["mean_rg"],
                    "min_rg_mmff": r["min_rg"],
                    "max_rg_mmff": r["max_rg"],
                    "energy_rg_correlation": r.get("energy_rg_corr"),
                }
            )

    with open("experiments/iter_25_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)

    print(final_text)
    print("\nFiles written:")
    print(" - experiments/iter_25_output.txt")
    print(" - experiments/iter_25_descriptors.json")


if __name__ == "__main__":
    main()
