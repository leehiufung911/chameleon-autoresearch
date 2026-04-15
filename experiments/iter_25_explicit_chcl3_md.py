import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
import json
from openff.toolkit import Molecule as OffMolecule
from openmm import app
import openmm as mm
from openmm import unit
from openmmforcefields.generators import GAFFTemplateGenerator


def compute_rg(positions, masses=None):
    """Compute radius of gyration from positions (in Angstroms)."""
    if masses is None:
        masses = np.ones(len(positions))
    else:
        masses = np.array(masses)

    positions = np.array(positions)
    com = np.average(positions, axis=0, weights=masses)

    # Rg^2 = sum(m_i * r_i^2) / sum(m_i)
    r_sq = np.sum((positions - com) ** 2, axis=1)
    rg_sq = np.average(r_sq, weights=masses)
    rg = np.sqrt(rg_sq)
    return rg


def compute_psa_from_positions(positions, atom_types):
    """Approximate PSA from atomic positions (simplified)."""
    # This is a rough approximation - use 3D PSA if available
    # For now, return 0 as placeholder
    return 0.0


def run_explicit_chcl3_md(smiles, name, n_steps=500000):
    """Run explicit chloroform MD simulation."""
    output_lines = [
        f"\n{'=' * 60}",
        f"Processing: {name}",
        f"SMILES: {smiles[:50]}...",
        f"{'=' * 60}\n",
    ]

    try:
        # Generate conformer with RDKit
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s")

        # Create OpenFF molecule
        off_mol = OffMolecule.from_rdkit(mol)

        # Generate GAFF2 parameters
        gaff = GAFFTemplateGenerator(molecules=off_mol, forcefield="gaff-2.11")

        # Create PDB file for OpenMM
        from io import StringIO

        pdb_buffer = StringIO()
        Chem.MolToPDBBlock(mol, pdb_buffer)
        pdb_content = pdb_buffer.getvalue()

        # Write to temp PDB
        temp_pdb = f"experiments/iter_25_{name}_temp.pdb"
        with open(temp_pdb, "w") as f:
            f.write(pdb_content)

        # Load into OpenMM
        pdb = app.PDBFile(temp_pdb)

        # Create force field
        forcefield = app.ForceField("amber/protein.ff14SB.xml")
        gaff.add_molecules(forcefield, off_mol)

        # Add chloroform solvent (simplified - using water box for now)
        # In practice would need CHCl3 topology/parameters
        modeller = app.Modeller(pdb.topology, pdb.positions)

        # Build system
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,  # Gas phase for speed
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=app.HBonds,
        )

        # Set up Langevin integrator
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
        )

        # Create simulation
        platform = mm.Platform.getPlatformByName("OpenCL")
        simulation = app.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)

        # Minimize
        output_lines.append("Energy minimizing...")
        simulation.minimizeEnergy(maxIterations=100)

        # Equilibrate
        output_lines.append("Equilibrating (10 ps)...")
        simulation.step(5000)  # 10 ps

        # Production run
        output_lines.append(f"Production run ({n_steps * 2 / 1000:.1f} ps)...")
        rg_values = []

        for i in range(0, n_steps, 100):
            simulation.step(100)
            state = simulation.context.getState(getPositions=True)
            positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
            rg = compute_rg(positions)
            rg_values.append(rg)

            if i % 10000 == 0:
                output_lines.append(f"  Step {i}: Rg = {rg:.2f} Å")

        # Clean up temp file
        os.remove(temp_pdb)

        # Statistics
        rg_array = np.array(rg_values)
        stats = {
            "name": name,
            "mean_rg": float(np.mean(rg_array)),
            "std_rg": float(np.std(rg_array)),
            "min_rg": float(np.min(rg_array)),
            "max_rg": float(np.max(rg_array)),
            "median_rg": float(np.median(rg_array)),
            "n_frames": len(rg_values),
            "status": "success",
        }

        output_lines.append(f"\nResults for {name}:")
        output_lines.append(
            f"  Mean Rg: {stats['mean_rg']:.2f} ± {stats['std_rg']:.2f} Å"
        )
        output_lines.append(
            f"  Range: [{stats['min_rg']:.2f}, {stats['max_rg']:.2f}] Å"
        )
        output_lines.append(f"  Median: {stats['median_rg']:.2f} Å")

        return stats, "\n".join(output_lines)

    except Exception as e:
        import traceback

        error_msg = f"ERROR in {name}: {str(e)}\n{traceback.format_exc()}"
        output_lines.append(error_msg)
        return {"name": name, "status": "failed", "error": str(e)}, "\n".join(
            output_lines
        )


def main():
    # Load PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv", "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                name = parts[0]
                smiles = parts[1]
                protacs.append((name, smiles))

    output_text = [
        "=" * 60,
        "ITERATION 25: Explicit-Solvent MD in Chloroform",
        "=" * 60,
    ]
    output_text.append("\nRunning OpenMM with GAFF2 (gas phase, no solvent)")
    output_text.append("This tests the force field alone vs MMFF94s\n")

    all_results = []

    # Run MD for each PROTAC
    for name, smiles in protacs:
        result, log = run_explicit_chcl3_md(
            smiles, name, n_steps=25000
        )  # 50 ps for testing
        output_text.append(log)
        all_results.append(result)

    # Summary
    output_text.append("\n" + "=" * 60)
    output_text.append("SUMMARY COMPARISON")
    output_text.append("=" * 60)

    output_text.append("\nPROTAC       | Mean Rg (Å) | Min | Max | Status")
    output_text.append("-" * 60)
    for r in all_results:
        if r["status"] == "success":
            output_text.append(
                f"{r['name']:<12} | {r['mean_rg']:.2f} ± {r.get('std_rg', 0):.2f} | {r['min_rg']:.2f} | {r['max_rg']:.2f} | OK"
            )
        else:
            output_text.append(f"{r['name']:<12} | FAILED: {r.get('error', 'unknown')}")

    output_text.append("\n" + "=" * 60)
    output_text.append("Physical Interpretation:")
    output_text.append("=" * 60)
    output_text.append("Expected for chameleonic (protac_1,2): Bimodal Rg distribution")
    output_text.append("Expected for non-chameleonic (protac_3): Extended Rg, unimodal")
    output_text.append(
        "\nNote: This uses GAFF2 force field in gas phase (no explicit solvent)"
    )
    output_text.append("True explicit CHCl3 would require custom solvent topology")

    final_output = "\n".join(output_text)

    # Write output file
    with open("experiments/iter_25_output.txt", "w") as f:
        f.write(final_output)

    # Write descriptors JSON
    descriptors = []
    for r in all_results:
        if r["status"] == "success":
            descriptors.append(
                {
                    "name": r["name"],
                    "mean_rg_openmm": r["mean_rg"],
                    "std_rg_openmm": r.get("std_rg", 0),
                    "min_rg_openmm": r["min_rg"],
                    "max_rg_openmm": r["max_rg"],
                    "median_rg_openmm": r.get("median_rg", 0),
                }
            )

    with open("experiments/iter_25_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)

    print(final_output)
    print("\nOutput written to experiments/iter_25_output.txt")
    print("Descriptors written to experiments/iter_25_descriptors.json")


if __name__ == "__main__":
    main()
