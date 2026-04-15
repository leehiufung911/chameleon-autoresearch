import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys

sys.path.insert(0, "chameleon_local")

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import json

# Position-Based Dynamics (PBD) Implementation
# Based on NVIDIA PhysX / game physics constraints


def compute_rg(positions, masses):
    """Compute radius of gyration from positions and masses"""
    center = np.average(positions, axis=0, weights=masses)
    deviations = positions - center
    rg_squared = np.average(np.sum(deviations**2, axis=1), weights=masses)
    return np.sqrt(rg_squared)


class PBDSolver:
    """Minimal Position-Based Dynamics solver for molecular conformers"""

    def __init__(self, mol, positions, stiffness=0.5, iterations=50):
        self.mol = mol
        self.n_atoms = mol.GetNumAtoms()
        self.stiffness = stiffness
        self.iterations = iterations
        self.constraints = []

        # Initial positions
        self.positions = positions.astype(np.float64).copy()
        self.masses = np.array(
            [max(1.0, atom.GetMass()) for atom in mol.GetAtoms()], dtype=np.float64
        )

    def add_minimum_distance_constraint(self, atom_i, atom_j, min_distance):
        """Add minimum distance constraint (only enforced if violated)"""
        self.constraints.append(("min_distance", atom_i, atom_j, min_distance))

    def solve_min_distance(self, i, j, min_dist):
        """Project minimum distance constraint (push apart if too close)"""
        delta = self.positions[j] - self.positions[i]
        current_length = np.linalg.norm(delta)

        if current_length < min_dist and current_length > 1e-6:
            error = min_dist - current_length
            correction = delta / current_length * error * self.stiffness

            w_i = 1.0 / self.masses[i]
            w_j = 1.0 / self.masses[j]
            w_sum = w_i + w_j

            self.positions[i] -= correction * w_i / w_sum
            self.positions[j] += correction * w_j / w_sum

    def compute_rg(self):
        """Compute radius of gyration"""
        return compute_rg(self.positions, self.masses)

    def step(self):
        """One PBD iteration - solve all constraints"""
        for constraint in self.constraints:
            if constraint[0] == "min_distance":
                _, i, j, min_dist = constraint
                self.solve_min_distance(i, j, min_dist)

    def solve(self):
        """Run multiple PBD iterations"""
        for _ in range(self.iterations):
            self.step()
        return self.positions.copy()


def run_pbd_experiment(smiles, name, n_conformers=20):
    """Run PBD on a molecule and return metrics"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    rgs_mmff = []
    rgs_pbd = []

    for seed in range(n_conformers):
        # Generate initial conformer with ETKDG
        params = AllChem.ETKDGv3()
        params.randomSeed = seed + 12345
        params.pruneRmsThresh = 0.5

        mol_copy = Chem.Mol(mol)
        AllChem.EmbedMolecule(mol_copy, params)
        if mol_copy.GetNumConformers() == 0:
            continue

        # Get initial positions
        conf = mol_copy.GetConformer()
        positions = conf.GetPositions().astype(np.float64)
        masses = np.array(
            [max(1.0, atom.GetMass()) for atom in mol_copy.GetAtoms()], dtype=np.float64
        )

        # MMFF minimization
        AllChem.MMFFOptimizeMolecule(mol_copy, mmffVariant="MMFF94s", maxIters=500)

        # Measure MMFF state
        positions_mmff = mol_copy.GetConformer().GetPositions().astype(np.float64)
        rg_mmff = compute_rg(positions_mmff, masses)
        rgs_mmff.append(rg_mmff)

        # PBD: Apply minimum distance constraints to preserve extended states
        # Strategy: Add anti-collapse constraints that keep the molecule from
        # collapsing below a certain size. This simulates entropic resistance.

        solver = PBDSolver(mol_copy, positions_mmff, stiffness=0.2, iterations=20)

        # Find atoms in "warhead" regions (first and last heavy atoms in SMILES)
        n_atoms = mol_copy.GetNumAtoms()

        # Simple heuristic: constrain first 5 and last 5 atoms to maintain distance
        # This simulates the linker keeping warheads apart in alkyl-linked PROTACs
        n_end = min(5, n_atoms // 5)

        # Compute current end-to-end distance
        end_to_end = np.linalg.norm(positions_mmff[0] - positions_mmff[-1])

        # Add minimum distance constraints between warhead regions
        # This is the key PBD insight: enforce entropic extension
        for i in range(n_end):
            for j in range(n_atoms - n_end, n_atoms):
                dist = np.linalg.norm(positions_mmff[i] - positions_mmff[j])
                # Only constrain if already somewhat extended (> 10 A)
                if dist > 8.0:
                    # Maintain at least 70% of current distance
                    solver.add_minimum_distance_constraint(i, j, dist * 0.7)

        # Solve constraints
        solver.solve()

        # Measure PBD state
        rg_pbd = solver.compute_rg()
        rgs_pbd.append(rg_pbd)

    if len(rgs_mmff) == 0:
        return None

    return {
        "name": name,
        "rg_mmff_mean": float(np.mean(rgs_mmff)),
        "rg_mmff_std": float(np.std(rgs_mmff)),
        "rg_mmff_max": float(np.max(rgs_mmff)),
        "rg_pbd_mean": float(np.mean(rgs_pbd)) if rgs_pbd else 0.0,
        "rg_pbd_std": float(np.std(rgs_pbd)) if rgs_pbd else 0.0,
        "rg_pbd_max": float(np.max(rgs_pbd)) if rgs_pbd else 0.0,
        "rg_extension_ratio": float(np.mean(rgs_pbd) / np.mean(rgs_mmff))
        if rgs_mmff and rgs_pbd
        else 0.0,
        "rg_diff": float(np.mean(rgs_pbd) - np.mean(rgs_mmff))
        if rgs_mmff and rgs_pbd
        else 0.0,
        "n_conformers": len(rgs_mmff),
    }


def main():
    # Read user PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                protacs.append((parts[0], parts[1]))

    results = []
    output_lines = []

    output_lines.append("=" * 70)
    output_lines.append(
        "ITERATION 17: Position-Based Dynamics (PBD) for Extended Conformers"
    )
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append(
        "Hypothesis: PBD anti-collapse constraints preserve extended conformations"
    )
    output_lines.append(
        "that MMFF94s collapses, simulating entropic resistance to folding in"
    )
    output_lines.append("alkyl-linked PROTACs like protac_3.")
    output_lines.append("")

    for name, smiles in protacs:
        output_lines.append(f"\nProcessing {name}...")
        result = run_pbd_experiment(smiles, name, n_conformers=15)
        if result:
            results.append(result)
            output_lines.append(f"  N conformers: {result['n_conformers']}")
            output_lines.append(
                f"  MMFF Rg:  {result['rg_mmff_mean']:.2f} +/- {result['rg_mmff_std']:.2f} A (max: {result['rg_mmff_max']:.2f})"
            )
            output_lines.append(
                f"  PBD Rg:   {result['rg_pbd_mean']:.2f} +/- {result['rg_pbd_std']:.2f} A (max: {result['rg_pbd_max']:.2f})"
            )
            output_lines.append(
                f"  Extension ratio: {result['rg_extension_ratio']:.3f}"
            )
            output_lines.append(f"  Rg difference: {result['rg_diff']:+.2f} A")
        else:
            output_lines.append(f"  FAILED to process {name}")

    # Summary
    output_lines.append("\n" + "=" * 70)
    output_lines.append("SUMMARY")
    output_lines.append("=" * 70)

    for r in results:
        output_lines.append(f"\n{r['name']}:")
        output_lines.append(f"  Rg extension via PBD: {r['rg_diff']:+.2f} A")
        output_lines.append(f"  Ratio (PBD/MMFF): {r['rg_extension_ratio']:.3f}")

    # Physical interpretation
    output_lines.append("\n" + "=" * 70)
    output_lines.append("PHYSICAL INTERPRETATION")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append(
        "The PBD anti-collapse constraint simulates entropic resistance to folding."
    )
    output_lines.append(
        "In explicit CHCl3 solvent, alkyl linkers (protac_3) stay extended due to"
    )
    output_lines.append(
        "conformational entropy. In implicit solvent MMFF, they collapse."
    )
    output_lines.append("")
    output_lines.append("Prediction for chameleonic molecules (protac_1, protac_2):")
    output_lines.append("  - Small Rg change (already compact in both methods)")
    output_lines.append("  - Ratio near 1.0 or slightly higher")
    output_lines.append("")
    output_lines.append("Prediction for non-chameleonic alkyl-linked (protac_3):")
    output_lines.append("  - Larger Rg increase (PBD corrects MMFF collapse)")
    output_lines.append("  - Ratio > 1.0 (PBD maintains extended conformers)")
    output_lines.append("")
    output_lines.append(
        "Key insight: If PBD makes protac_3 more extended while keeping"
    )
    output_lines.append("protac_1/2 similar, then constraint-based physics is a valid")
    output_lines.append("alternative to full explicit-solvent MD.")

    output = "\n".join(output_lines)

    # Write descriptors for all 49 molecules (simplified - just PROTACs for now)
    descriptors = []
    for r in results:
        descriptors.append(
            {
                "name": r["name"],
                "rg_mmff": r["rg_mmff_mean"],
                "rg_pbd": r["rg_pbd_mean"],
                "rg_ratio": r["rg_extension_ratio"],
                "rg_diff": r["rg_diff"],
            }
        )

    with open("experiments/iter_17_descriptors.json", "w") as f:
        json.dump(descriptors, f, indent=2)

    with open("experiments/iter_17_output.txt", "w") as f:
        f.write(output)
    print(output, flush=True)


if __name__ == "__main__":
    main()
