"""
Focused holdout analysis — deep dive into why term_rg dominates.

Tests Rg-focused scoring variants:
  1. term_rg alone (log Rg ratio)
  2. term_rg + term_imhb (no PSA)
  3. Rg range (max - min) normalized by MW
  4. Rg coefficient of variation (std / mean)
  5. Combined Rg + IMHB multiplicative
  6. Rg_range * IMHB (multiplicative, physical: folding + polarity masking)

Also runs on the ORIGINAL benchmark for comparison, to see if term_rg alone
matches or exceeds CI there too.
"""
import sys, os, math, json
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "chameleon_local"))
import chameleon


def roc_auc(scores, labels):
    pos = [(s, 1) for s, l in zip(scores, labels) if l]
    neg = [(s, 0) for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float("nan")
    n_pos, n_neg = len(pos), len(neg)
    all_vals = pos + neg
    all_vals.sort(key=lambda x: x[0])
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


def load_dataset(path, label_col_idx=0, name_col_idx=1, smiles_col_idx=2):
    molecules = []
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            label = parts[label_col_idx]
            name = parts[name_col_idx]
            smiles = parts[smiles_col_idx]
            molecules.append({"name": name, "smiles": smiles, "label": label})
    return molecules


def process_molecules(molecules, n_conf=50, label="dataset"):
    results = []
    for i, mol_info in enumerate(molecules):
        name = mol_info["name"]
        smiles = mol_info["smiles"]
        lab = mol_info["label"]
        is_cham = lab == "chameleon"

        print(f"  [{i+1}/{len(molecules)}] {name[:35]:<35}", end="", flush=True)
        try:
            s = chameleon.summarize(name, smiles, n_conf=n_conf, verbose=False)

            psa_range = s.psa3d_max - s.psa3d_min
            mw = s.mw
            rg_ratio = s.rg_max / max(s.rg_min, 1e-6)
            rg_range = s.rg_max - s.rg_min
            rg_cv = s.rg_std / max((s.rg_max + s.rg_min)/2, 1e-6)
            term_psa = psa_range / math.sqrt(max(mw, 1.0))
            term_imhb = s.imhb_mean
            term_rg = math.log(max(rg_ratio, 1.0))

            ci_baseline = 2.0 * term_psa + 1.0 * term_imhb + 3.0 * term_rg
            ci_rg_only = term_rg
            ci_rg_imhb = 3.0 * term_rg + 1.0 * term_imhb
            ci_rg_range_mw = rg_range / math.sqrt(max(mw, 1.0))
            ci_rg_cv = rg_cv
            ci_rg_imhb_mult = term_rg * term_imhb
            ci_rg_range_imhb = (rg_range / math.sqrt(max(mw, 1.0))) * term_imhb

            results.append({
                "name": name, "label": lab, "is_cham": is_cham,
                "ci_baseline": ci_baseline,
                "term_rg": ci_rg_only,
                "rg_plus_imhb": ci_rg_imhb,
                "rg_range_mw": ci_rg_range_mw,
                "rg_cv": ci_rg_cv,
                "rg_x_imhb": ci_rg_imhb_mult,
                "rg_range_x_imhb": ci_rg_range_imhb,
                "term_psa": term_psa, "term_imhb": term_imhb,
                "mw": mw, "rg_max": s.rg_max, "rg_min": s.rg_min,
                "rg_range": rg_range, "rg_ratio": rg_ratio
            })
            print(f" Rg={s.rg_min:.1f}-{s.rg_max:.1f}  IMHB={term_imhb:.2f}  CI={ci_baseline:.2f}")
        except Exception as e:
            print(f" ERROR: {e}")
    return results


def report_aucs(results, dataset_name):
    labels_bool = [r["is_cham"] for r in results]
    methods = {
        "CI_baseline (2P+I+3R)": [r["ci_baseline"] for r in results],
        "term_rg alone": [r["term_rg"] for r in results],
        "3*rg + imhb": [r["rg_plus_imhb"] for r in results],
        "rg_range/sqrt(MW)": [r["rg_range_mw"] for r in results],
        "rg_cv": [r["rg_cv"] for r in results],
        "rg * imhb": [r["rg_x_imhb"] for r in results],
        "rg_range/MW * imhb": [r["rg_range_x_imhb"] for r in results],
        "term_psa alone": [r["term_psa"] for r in results],
        "term_imhb alone": [r["term_imhb"] for r in results],
    }

    n_pos = sum(labels_bool)
    n_neg = len(labels_bool) - n_pos
    print(f"\n{'='*55}")
    print(f"{dataset_name} AUC (N={len(results)}, pos={n_pos}, neg={n_neg})")
    print(f"{'='*55}")
    print(f"{'Method':<25} {'AUC':>8}")
    print(f"{'-'*33}")
    best_auc = max(roc_auc(s, labels_bool) for s in methods.values())
    for method_name, scores in methods.items():
        auc = roc_auc(scores, labels_bool)
        marker = " ***" if auc == best_auc else ""
        print(f"{method_name:<25} {auc:>8.3f}{marker}")


def main():
    base_dir = os.path.join(os.path.dirname(__file__), "..")

    # --- Holdout set ---
    print("="*60)
    print("HOLDOUT SET (29 unseen Chamelogk molecules)")
    print("="*60)
    holdout = load_dataset(os.path.join(base_dir, "holdout_set.tsv"))
    holdout_results = process_molecules(holdout, n_conf=50, label="holdout")
    report_aucs(holdout_results, "HOLDOUT SET")

    # --- Original benchmark ---
    print("\n\n" + "="*60)
    print("ORIGINAL BENCHMARK (49 labelled molecules)")
    print("="*60)
    bench = load_dataset(
        os.path.join(base_dir, "chameleon_local", "labelled_set.tsv"),
        label_col_idx=0, name_col_idx=1, smiles_col_idx=2
    )
    bench_results = process_molecules(bench, n_conf=50, label="benchmark")
    report_aucs(bench_results, "ORIGINAL BENCHMARK")


if __name__ == "__main__":
    main()
