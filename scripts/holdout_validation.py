"""
Holdout validation — test promising hypotheses from the research loop on 29
unseen bRo5 compounds from the Chamelogk paper (never shown to Gemini).

Tests:
  1. Baseline CI (2*term_psa + 1*term_imhb + 3*term_rg)
  2. Multiplicative CI coupling (term_psa + 3*term_rg) * term_imhb
  3. MaxCPath >= 5 & IMHB < 1.0 veto (halves CI for aliphatic-linker molecules)
  4. Combined: multiplicative CI + MaxCPath veto
  5. H-Bond Saturation Gain (iter 9): IMHB / HBD as normalization

Reports AUC for each method on the holdout set.
"""
import sys, os, math, json
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "chameleon_local"))
import chameleon

# ---------------------------------------------------------------------------
# Helper: longest contiguous non-ring carbon path (from iter 13)
# ---------------------------------------------------------------------------
def get_max_nonring_carbon_path(mol):
    nonring_c = [a.GetIdx() for a in mol.GetAtoms()
                 if a.GetSymbol() == "C" and not a.IsInRing()]
    if not nonring_c:
        return 0
    n = len(nonring_c)
    idx_map = {idx: i for i, idx in enumerate(nonring_c)}
    adj = [[False]*n for _ in range(n)]
    for i, idx in enumerate(nonring_c):
        atom = mol.GetAtomWithIdx(idx)
        for nb in atom.GetNeighbors():
            if nb.GetIdx() in idx_map:
                adj[i][idx_map[nb.GetIdx()]] = True
    # BFS from each node to find diameter
    from collections import deque
    max_path = 0
    for start in range(n):
        visited = [False]*n
        visited[start] = True
        q = deque([(start, 0)])
        while q:
            node, dist = q.popleft()
            max_path = max(max_path, dist)
            for j in range(n):
                if adj[node][j] and not visited[j]:
                    visited[j] = True
                    q.append((j, dist+1))
    return max_path + 1  # path length = edges + 1


def count_hbd(mol):
    """Count H-bond donors (NH, OH)."""
    return Descriptors.NumHDonors(mol)


def roc_auc(scores, labels):
    """Mann-Whitney U AUC (no sklearn dependency)."""
    pos = [(s, 1) for s, l in zip(scores, labels) if l]
    neg = [(s, 0) for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float("nan")
    n_pos, n_neg = len(pos), len(neg)
    all_vals = pos + neg
    all_vals.sort(key=lambda x: x[0])
    # assign ranks
    ranks = [0.0] * len(all_vals)
    i = 0
    while i < len(all_vals):
        j = i
        while j + 1 < len(all_vals) and all_vals[j+1][0] == all_vals[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j+1):
            ranks[k] = avg_rank
        i = j + 1
    sum_pos_ranks = sum(r for r, (_, lab) in zip(ranks, all_vals) if lab == 1)
    u = sum_pos_ranks - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


def main():
    # Load holdout set
    holdout_path = os.path.join(os.path.dirname(__file__), "..", "holdout_set.tsv")
    molecules = []
    with open(holdout_path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            label, name, smiles = parts[0], parts[1], parts[2]
            molecules.append({"name": name, "smiles": smiles, "label": label})

    print(f"Loaded {len(molecules)} holdout molecules "
          f"({sum(1 for m in molecules if m['label']=='chameleon')} cham, "
          f"{sum(1 for m in molecules if m['label']=='nonchameleon')} non-cham)")
    print()

    N_CONF = 50  # reasonable for validation
    results = []

    for i, mol_info in enumerate(molecules):
        name = mol_info["name"]
        smiles = mol_info["smiles"]
        label = mol_info["label"]
        is_cham = label == "chameleon"

        print(f"[{i+1}/{len(molecules)}] {name[:30]:<30} ({label})", end="", flush=True)
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(" -> SKIP (invalid SMILES)")
                continue

            s = chameleon.summarize(name, smiles, n_conf=N_CONF, verbose=False)

            # Standard CI components
            psa_range = s.psa3d_max - s.psa3d_min
            mw = s.mw
            rg_ratio = s.rg_max / max(s.rg_min, 1e-6)
            term_psa = psa_range / math.sqrt(max(mw, 1.0))
            term_imhb = s.imhb_mean
            term_rg = math.log(max(rg_ratio, 1.0))

            # 1. Baseline CI
            ci_baseline = 2.0 * term_psa + 1.0 * term_imhb + 3.0 * term_rg

            # 2. Multiplicative CI
            ci_mult = (term_psa + 3.0 * term_rg) * term_imhb

            # 3. MaxCPath veto on baseline CI
            max_c = get_max_nonring_carbon_path(mol)
            ci_veto = ci_baseline
            if max_c >= 5 and term_imhb < 1.0:
                ci_veto *= 0.5

            # 4. Combined: multiplicative + veto
            ci_combined = ci_mult
            if max_c >= 5 and term_imhb < 1.0:
                ci_combined *= 0.5

            # 5. H-Bond Saturation Gain (iter 9 idea)
            hbd = count_hbd(mol)
            sat_gain = term_imhb / max(hbd, 1) if hbd > 0 else 0.0
            ci_sat = 2.0 * term_psa + 2.0 * sat_gain + 3.0 * term_rg

            results.append({
                "name": name, "label": label, "is_cham": is_cham,
                "ci_baseline": ci_baseline,
                "ci_mult": ci_mult,
                "ci_veto": ci_veto,
                "ci_combined": ci_combined,
                "ci_sat": ci_sat,
                "term_psa": term_psa, "term_imhb": term_imhb,
                "term_rg": term_rg, "max_c": max_c, "hbd": hbd,
                "mw": mw, "psa_range": psa_range
            })
            print(f"  CI={ci_baseline:.2f}  CI_mult={ci_mult:.2f}  "
                  f"MaxC={max_c}  IMHB={term_imhb:.2f}")

        except Exception as e:
            print(f" -> ERROR: {e}")
            continue

    if not results:
        print("No results computed!")
        return

    # --- Compute AUCs ---
    labels_bool = [r["is_cham"] for r in results]
    methods = {
        "CI_baseline": [r["ci_baseline"] for r in results],
        "CI_multiplicative": [r["ci_mult"] for r in results],
        "CI_maxcpath_veto": [r["ci_veto"] for r in results],
        "CI_combined": [r["ci_combined"] for r in results],
        "CI_saturation": [r["ci_sat"] for r in results],
        # Also test individual channels
        "term_psa": [r["term_psa"] for r in results],
        "term_imhb": [r["term_imhb"] for r in results],
        "term_rg": [r["term_rg"] for r in results],
        "psa_range": [r["psa_range"] for r in results],
    }

    n_pos = sum(labels_bool)
    n_neg = len(labels_bool) - n_pos
    print(f"\n{'='*70}")
    print(f"HOLDOUT SET AUC RESULTS (N={len(results)}, pos={n_pos}, neg={n_neg})")
    print(f"{'='*70}")
    print(f"{'Method':<25} {'AUC':>8}")
    print(f"{'-'*33}")
    for method_name, scores in methods.items():
        auc = roc_auc(scores, labels_bool)
        marker = " <--- BEST" if auc == max(roc_auc(s, labels_bool) for s in methods.values()) else ""
        print(f"{method_name:<25} {auc:>8.3f}{marker}")

    # --- Per-molecule detail ---
    print(f"\n{'='*70}")
    print(f"PER-MOLECULE DETAIL")
    print(f"{'='*70}")
    print(f"{'Name':<30} {'Label':<8} {'CI_base':>7} {'CI_mult':>7} "
          f"{'CI_veto':>7} {'CI_comb':>7} {'MaxC':>4} {'IMHB':>5}")
    print("-" * 85)
    for r in sorted(results, key=lambda x: -x["ci_baseline"]):
        print(f"{r['name'][:30]:<30} "
              f"{'CHAM' if r['is_cham'] else 'NON':>5}   "
              f"{r['ci_baseline']:>7.2f} {r['ci_mult']:>7.2f} "
              f"{r['ci_veto']:>7.2f} {r['ci_combined']:>7.2f} "
              f"{r['max_c']:>4} {r['term_imhb']:>5.2f}")

    # Save results JSON
    out_path = os.path.join(os.path.dirname(__file__), "..", "holdout_validation_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved results to {out_path}")


if __name__ == "__main__":
    main()
