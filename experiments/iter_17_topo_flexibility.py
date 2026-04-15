import os

os.environ["PYTHONUNBUFFERED"] = "1"

import sys

sys.path.insert(0, "chameleon_local")

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
import json

# Topological Flexibility Analysis
# Inspired by rigidity theory and polymer physics
# Question: Can we predict chameleonicity from 2D topology alone?


def count_ether_oxygens(mol):
    """Count ether oxygen atoms (excluding carbonyls)"""
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O":
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            # Ether: O connected to 2 carbons (not carbonyl)
            if neighbors == ["C", "C"]:
                count += 1
    return count


def count_alkyl_chain_carbons(mol):
    """Count carbons in long alkyl chains (>3 consecutive CH2)"""
    count = 0
    visited = set()

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetIdx() not in visited:
            # Check for chain of sp3 carbons
            if atom.GetHybridization() == Chem.HybridizationType.SP3:
                chain = []
                current = atom
                while current and current.GetSymbol() == "C":
                    if current.GetIdx() in visited:
                        break
                    visited.add(current.GetIdx())
                    chain.append(current)
                    # Find next sp3 carbon
                    next_c = None
                    for n in current.GetNeighbors():
                        if (
                            n.GetSymbol() == "C"
                            and n.GetHybridization() == Chem.HybridizationType.SP3
                        ):
                            if n.GetIdx() not in visited:
                                next_c = n
                                break
                    current = next_c

                if len(chain) > 3:
                    count += len(chain)

    return count


def find_shortest_linker_path(mol):
    """
    Find the shortest path between the two largest rings or between
    the two most distant functional groups
    """
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    if len(rings) >= 2:
        # Find path between two largest rings
        ring1 = max(rings, key=len)
        ring2 = [r for r in rings if r != ring1]
        if ring2:
            ring2 = max(ring2, key=len)
            # Find shortest path between rings
            # Simplified: count atoms between them
            return estimate_linker_length(mol, ring1, ring2)

    return 0


def estimate_linker_length(mol, atoms1, atoms2):
    """Estimate linker length between two sets of atoms"""
    # Find minimum distance between any atom in set1 and set2
    min_dist = float("inf")
    for i in atoms1:
        for j in atoms2:
            try:
                path = Chem.rdmolops.GetShortestPath(mol, i, j)
                if path:
                    dist = len(path) - 1  # Number of bonds
                    if dist < min_dist and dist > 0:
                        min_dist = dist
            except:
                pass

    return min_dist if min_dist != float("inf") else 0


def compute_flexibility_entropy(mol):
    """
    Compute a 2D flexibility score based on rotatable bonds.

    Key insight: Chameleonicity requires enough flexibility to fold,
    but not so much that the molecule is just floppy.

    Formula:
    - Count rotatable bonds
    - Penalize very long alkyl chains (high conformational entropy = hard to fold)
    - Reward PEG linkers (ether oxygens = gauche preference = easier folding)
    """
    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_ether = count_ether_oxygens(mol)
    n_alkyl = count_alkyl_chain_carbons(mol)
    mw = Descriptors.MolWt(mol)

    # Topological flexibility index
    # High rotatable bonds = flexible
    # Ether oxygens = promote folding (gauche effect)
    # Long alkyl chains = entropic penalty (anti preference)

    if n_rot == 0:
        return 0.0, {}

    # Ether density in rotatable bonds
    ether_ratio = n_ether / n_rot if n_rot > 0 else 0

    # Alkyl penalty (long chains resist folding)
    alkyl_penalty = min(1.0, n_alkyl / 10.0)

    # Flexibility-entropy balance
    # High flexibility + ether-rich = chameleonic
    # High flexibility + alkyl-rich = just floppy (non-chameleonic)

    flexibility = n_rot / np.sqrt(mw / 100)  # Size-normalized flexibility

    # Entropic folding cost: high for alkyl, low for ether
    entropic_cost = flexibility * (1 - ether_ratio) * (1 + alkyl_penalty)

    # Chameleonic potential: flexibility minus entropic cost
    cham_potential = flexibility - entropic_cost

    details = {
        "n_rotatable": n_rot,
        "n_ether_o": n_ether,
        "n_alkyl_c": n_alkyl,
        "mw": mw,
        "ether_ratio": ether_ratio,
        "alkyl_penalty": alkyl_penalty,
        "flexibility": flexibility,
        "entropic_cost": entropic_cost,
        "cham_potential": cham_potential,
    }

    return cham_potential, details


def compute_graph_rigidity(mol):
    """
    Use 2D graph theory concepts from rigidity theory.

    Pebble game insight: A graph is rigid if it has enough
    constraints. Molecules with rings are locally rigid.

    For chameleonicity: We want the LINKER between rigid parts
    to be flexible but not too flexible.
    """
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    n_rings = len(rings)

    if n_rings < 2:
        return 0.0, {"n_rings": n_rings, "linker_bonds": 0}

    # Find shortest path between two largest rings
    ring_sizes = [(len(r), r) for r in rings]
    ring_sizes.sort(reverse=True)

    if len(ring_sizes) >= 2:
        ring1 = ring_sizes[0][1]
        ring2 = ring_sizes[1][1]

        # Find shortest path
        min_path_len = float("inf")
        for i in ring1:
            for j in ring2:
                try:
                    path = Chem.rdmolops.GetShortestPath(mol, i, j)
                    if path:
                        path_len = (
                            len(path)
                            - len(set(path) & set(ring1))
                            - len(set(path) & set(ring2))
                        )
                        if path_len < min_path_len and path_len > 0:
                            min_path_len = path_len
                except:
                    pass

        if min_path_len == float("inf"):
            min_path_len = 0
    else:
        min_path_len = 0

    # Rigidity score
    # Many rings = rigid scaffold (good for chameleonicity)
    # Long linker = flexible (good) but too long = floppy (bad)

    ring_rigidity = min(5.0, n_rings)
    linker_flex = min(3.0, min_path_len / 5.0) if min_path_len > 0 else 0

    # Optimal: rigid rings + medium linker
    rigidity_score = ring_rigidity + linker_flex - abs(linker_flex - 1.5)

    return rigidity_score, {
        "n_rings": n_rings,
        "linker_bonds": min_path_len,
        "ring_rigidity": ring_rigidity,
        "linker_flex": linker_flex,
    }


def analyze_molecule(smiles, name):
    """Full topological analysis"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Compute descriptors
    cham_pot, flex_details = compute_flexibility_entropy(mol)
    rigidity, rigid_details = compute_graph_rigidity(mol)

    # Combined topological chameleonicity index
    # Weight: 70% flexibility-entropy, 30% rigidity
    topo_ci = 0.7 * cham_pot + 0.3 * rigidity

    return {
        "name": name,
        "smiles": smiles,
        "topo_cham_potential": cham_pot,
        "topo_rigidity": rigidity,
        "topo_ci": topo_ci,
        **flex_details,
        **{f"rig_{k}": v for k, v in rigid_details.items()},
    }


def main():
    output_lines = []

    output_lines.append("=" * 70)
    output_lines.append("ITERATION 17: 2D Topological Flexibility Analysis")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append("Hypothesis: Chameleonicity can be predicted from 2D topology")
    output_lines.append("using rigidity theory and polymer physics principles, without")
    output_lines.append("any 3D conformer generation.")
    output_lines.append("")
    output_lines.append("Key insights from polymer physics:")
    output_lines.append(
        "  - PEG linkers: gauche conformations favor folding (low entropy cost)"
    )
    output_lines.append(
        "  - Alkyl linkers: anti conformations resist folding (high entropy cost)"
    )
    output_lines.append("  - Ring systems provide rigid scaffolds for IMHB formation")
    output_lines.append("  - Optimal chameleon: rigid rings + ether-rich medium linker")
    output_lines.append("")

    # Read user PROTACs
    protacs = []
    with open("chameleon_local/user_protacs.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                protacs.append((parts[0], parts[1]))

    results = []

    output_lines.append("-" * 70)
    output_lines.append("USER PROTACS ANALYSIS")
    output_lines.append("-" * 70)
    output_lines.append("")

    for name, smiles in protacs:
        result = analyze_molecule(smiles, name)
        if result:
            results.append(result)
            output_lines.append(f"{name}:")
            output_lines.append(f"  Rotatable bonds: {result['n_rotatable']}")
            output_lines.append(f"  Ether oxygens: {result['n_ether_o']}")
            output_lines.append(f"  Alkyl chain C: {result['n_alkyl_c']}")
            output_lines.append(f"  Ether ratio: {result['ether_ratio']:.3f}")
            output_lines.append(f"  Alkyl penalty: {result['alkyl_penalty']:.3f}")
            output_lines.append(f"  Flexibility: {result['flexibility']:.3f}")
            output_lines.append(f"  Entropic cost: {result['entropic_cost']:.3f}")
            output_lines.append(
                f"  Chameleonic potential: {result['topo_cham_potential']:.3f}"
            )
            output_lines.append(f"  Ring count: {result['rig_n_rings']}")
            output_lines.append(f"  Linker bonds: {result['rig_linker_bonds']}")
            output_lines.append(f"  TOPOLOGICAL CI: {result['topo_ci']:.3f}")
            output_lines.append("")

    # Analyze benchmark
    output_lines.append("-" * 70)
    output_lines.append("BENCHMARK MOLECULES (49 molecules)")
    output_lines.append("-" * 70)
    output_lines.append("")

    benchmark = []
    with open("chameleon_local/benchmark.tsv") as f:
        header = next(f).strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                name = parts[0]
                smiles = parts[1]
                label = parts[2]
                benchmark.append((name, smiles, label))

    benchmark_results = []
    for name, smiles, label in benchmark:
        result = analyze_molecule(smiles, name)
        if result:
            result["label"] = label
            benchmark_results.append(result)

    # Compute statistics by label
    cham_values = [
        r["topo_ci"] for r in benchmark_results if r["label"] == "chameleonic"
    ]
    noncham_values = [
        r["topo_ci"] for r in benchmark_results if r["label"] == "non_chameleonic"
    ]

    if cham_values and noncham_values:
        output_lines.append(
            f"Chameleonic (n={len(cham_values)}): mean={np.mean(cham_values):.3f}, std={np.std(cham_values):.3f}"
        )
        output_lines.append(
            f"Non-chameleonic (n={len(noncham_values)}): mean={np.mean(noncham_values):.3f}, std={np.std(noncham_values):.3f}"
        )
        output_lines.append(
            f"Separation: {np.mean(cham_values) - np.mean(noncham_values):.3f}"
        )

    output_lines.append("")

    # Qualitative check for PROTACs
    output_lines.append("-" * 70)
    output_lines.append("QUALITATIVE CHECK")
    output_lines.append("-" * 70)
    output_lines.append("")

    # Find CI values for PROTACs
    protac_cis = [(r["name"], r["topo_ci"]) for r in results]
    protac_cis.sort(key=lambda x: x[1], reverse=True)

    output_lines.append("PROTACs ranked by Topological CI:")
    for name, ci in protac_cis:
        gt = "CHAM" if name in ["protac_1", "protac_2"] else "NON"
        output_lines.append(f"  {name}: {ci:.3f} (ground truth: {gt})")

    # Check separation
    cham_scores = [ci for name, ci in protac_cis if name in ["protac_1", "protac_2"]]
    non_scores = [ci for name, ci in protac_cis if name == "protac_3"]

    if cham_scores and non_scores:
        output_lines.append("")
        output_lines.append(f"Chameleonic PROTACs avg CI: {np.mean(cham_scores):.3f}")
        output_lines.append(f"Non-chameleonic avg CI: {np.mean(non_scores):.3f}")
        output_lines.append(
            f"Difference: {np.mean(cham_scores) - np.mean(non_scores):+.3f}"
        )

        if np.mean(cham_scores) > np.mean(non_scores):
            output_lines.append("\n✓ Correct qualitative ordering!")
        else:
            output_lines.append("\n✗ Wrong ordering - method fails")

    output_lines.append("")
    output_lines.append("=" * 70)
    output_lines.append("PHYSICAL INTERPRETATION")
    output_lines.append("=" * 70)
    output_lines.append("")
    output_lines.append("The topological approach uses 2D insights:")
    output_lines.append(
        "1. Rotatable bonds = flexibility (required for chameleonicity)"
    )
    output_lines.append(
        "2. Ether oxygens = low entropic cost of folding (gauche effect)"
    )
    output_lines.append("3. Alkyl chains = high entropic cost (anti preference)")
    output_lines.append("4. Rings = rigid scaffolds for IMHB networks")
    output_lines.append("")
    output_lines.append("This is fundamentally different from 3D MMFF approaches.")
    output_lines.append(
        "It does not require conformer generation or energy minimization."
    )
    output_lines.append(
        "Trade-off: Less physical detail, but no implicit solvent bias."
    )

    output = "\n".join(output_lines)

    # Write output
    with open("experiments/iter_17_output.txt", "w") as f:
        f.write(output)
    print(output, flush=True)

    # Write descriptors
    with open("experiments/iter_17_descriptors.json", "w") as f:
        json.dump(benchmark_results, f, indent=2)

    print(
        f"\nWrote {len(benchmark_results)} molecule descriptors to iter_17_descriptors.json"
    )


if __name__ == "__main__":
    main()
