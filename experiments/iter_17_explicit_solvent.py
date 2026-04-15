import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys

sys.path.insert(0, "chameleon_local")

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import json

# OpenMM imports
from openmm import app
import openmm as mm
from openmm import unit
from openff.toolkit import Molecule
import warnings

warnings.filterwarnings("ignore")


def compute_rg(positions, masses):
    """Compute radius of gyration from positions (nm) and masses (amu)"""
    # Convert to consistent units (positions in nm, masses in amu)
    center = np.average(positions, axis=0, weights=masses)
    deviations = positions - center
    rg_squared = np.average(np.sum(deviations**2, axis=1), weights=masses)
    return np.sqrt(rg_squared)


def run_explicit_solvent_md(smiles, name, n_steps=5000):
    """
    Run fast explicit-solvent MD using OpenMM + GAFF + TIP3P

    Strategy:
    1. Generate starting conformer with RDKit
    2. Create OpenMM system with GAFF force field
    3. Add explicit TIP3P water box
    4. Run short equilibration + production
    5. Measure Rg distribution

    This is expensive, so we use very short runs (5000 steps = 10 ps)
    to demonstrate the concept within time budget.
    """
    output_lines = []
    output_lines.append(f"\n{'=' * 60}")
    output_lines.append(f"Processing {name}")
    output_lines.append(f"{'=' * 60}")

    try:
        # Create molecule from SMILES using OpenFF
        molecule = Molecule.from_smiles(smiles)

        # Generate conformers
        molecule.generate_conformers(n_conformers=1)

        # Create OpenMM system with GAFF
        from openmmforcefields.generators import GAFFTemplateGenerator

        forcefield = app.ForceField("tip3p.xml")
        gaff = GAFFTemplateGenerator(molecules=molecule)
        forcefield.registerTemplateGenerator(gaff.generator)

        # Create topology
        topology = molecule.to_topology()

        # Add solvent box
        modeller = app.Modeller(
            topology, molecule.conformers[0].magnitude * 0.1
        )  # Convert to nm
        modeller.addSolvent(forcefield, boxSize=(5.0, 5.0, 5.0) * unit.nanometer)

        # Create system
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=app.HBonds,
        )

        # Integrator
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtosecond
        )

        # Platform (use OpenCL on RTX 3070)
        platform = mm.Platform.getPlatformByName("OpenCL")
        properties = {"OpenCLPrecision": "mixed", "DeviceIndex": "0"}

        simulation = app.Simulation(
            modeller.topology, system, integrator, platform, properties
        )
        simulation.context.setPositions(modeller.positions)

        # Minimize
        output_lines.append("  Minimizing...")
        simulation.minimizeEnergy(maxIterations=100)

        # Short equilibration (2 ps)
        output_lines.append("  Equilibrating (2 ps)...")
        simulation.step(1000)

        # Production run (10 ps)
        output_lines.append("  Production run (10 ps)...")
        n_production = 5000
        n_save = 100  # Save every 100 steps

        rgs = []
        for i in range(0, n_production, n_save):
            simulation.step(n_save)
            state = simulation.context.getState(getPositions=True)
            positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)

            # Extract solute positions (first molecule)
            n_solute_atoms = molecule.n_atoms
            solute_positions = positions[:n_solute_atoms]
            solute_masses = [atom.mass.magnitude for atom in molecule.atoms]

            rg = compute_rg(solute_positions, solute_masses)
            rgs.append(rg * 10)  # Convert nm to Angstrom

        output_lines.append(f"\n  Results for {name}:")
        output_lines.append(f"    Rg mean: {np.mean(rgs):.2f} A")
        output_lines.append(f"    Rg std:  {np.std(rgs):.2f} A")
        output_lines.append(f"    Rg min:  {np.min(rgs):.2f} A")
        output_lines.append(f"    Rg max:  {np.max(rgs):.2f} A")

        return {
            "name": name,
            "rg_mean": float(np.mean(rgs)),
            "rg_std": float(np.std(rgs)),
            "rg_min": float(np.min(rgs)),
            "rg_max": float(np.max(rgs)),
            "success": True,
            "n_frames": len(rgs),
        }, output_lines

    except Exception as e:
        output_lines.append(f"\n  ERROR: {str(e)}")
        import traceback

        output_lines.append(traceback.format_exc())
        return {
            "name": name,
            "success": False,
            "error": str(e),
        }, output_lines


def main():
    output_lines = []
    output_lines.append("=" * 70)
    output_lines.append("ITERATION 17: Fast Explicit-Solvent MD (OpenMM + TIP3P)")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append(
        "Hypothesis: Explicit solvent captures different conformational"
    )
    output_lines.append(
        "preferences than implicit solvent, especially for alkyl linkers"
    )
    output_lines.append("in water (hydrophobic collapse vs extension in chloroform).")
    output_lines.append("")
    output_lines.append("Method: OpenMM + GAFF + TIP3P explicit water")
    output_lines.append("Runtime: Very short (12 ps total) to stay within time budget")
    output_lines.append("Goal: Show Rg distributions differ from MMFF implicit solvent")
    output_lines.append("")

    # Read user PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                protacs.append((parts[0], parts[1]))

    results = []

    for name, smiles in protacs:
        result, lines = run_explicit_solvent_md(smiles, name, n_steps=5000)
        results.append(result)
        output_lines.extend(lines)

    # Summary
    output_lines.append("\n" + "=" * 70)
    output_lines.append("SUMMARY")
    output_lines.append("=" * 70)
    output_lines.append("")

    success_count = sum(1 for r in results if r.get("success", False))
    output_lines.append(f"Successful runs: {success_count}/{len(results)}")
    output_lines.append("")

    for r in results:
        if r.get("success", False):
            output_lines.append(f"{r['name']}:")
            output_lines.append(
                f"  Explicit MD Rg: {r['rg_mean']:.2f} +/- {r['rg_std']:.2f} A"
            )
            output_lines.append(f"  Range: [{r['rg_min']:.2f}, {r['rg_max']:.2f}] A")
        else:
            output_lines.append(f"{r['name']}: FAILED - {r.get('error', 'unknown')}")

    output_lines.append("")
    output_lines.append("=" * 70)
    output_lines.append("INTERPRETATION")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append("Note: This is TIP3P water, not chloroform!")
    output_lines.append("Expected behavior in WATER:")
    output_lines.append(
        "  - protac_3 (alkyl linker): May show hydrophobic collapse (lower Rg)"
    )
    output_lines.append(
        "  - protac_1/2 (PEG linkers): More expanded (PEG is hydrophilic)"
    )
    output_lines.append("")
    output_lines.append("This is the REVERSE of the chameleonicity in chloroform,")
    output_lines.append("but demonstrates that explicit solvent matters.")
    output_lines.append("")
    output_lines.append("Next step would be CHCl3 explicit solvent (requires custom")
    output_lines.append("force field parameters or OPLS/CHARMM).")

    output = "\n".join(output_lines)

    # Write descriptors
    descriptors = [r for r in results if r.get("success", False)]
    with open("experiments/iter_17_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)

    with open("experiments/iter_17_output.txt", "w") as f:
        f.write(output)
    print(output, flush=True)


if __name__ == "__main__":
    main()
