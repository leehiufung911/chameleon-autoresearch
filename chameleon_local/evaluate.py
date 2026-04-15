"""
Merge a chameleon.py JSON results file with the literature-labelled set and
report classification performance for every scalar descriptor in the JSON.

Usage:
    py evaluate.py baseline_results.json labelled_set.tsv
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict


def read_labels(tsv_path: str) -> dict:
    labels = {}
    with open(tsv_path) as f:
        header = f.readline().strip().split("\t")
        iL = header.index("label")
        iN = header.index("name")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            labels[parts[iN]] = parts[iL]
    return labels


def roc_auc(scores, positives):
    """Mann-Whitney-U / rank-sum AUC (no sklearn).
    scores   : list[float]
    positives: list[bool]  -- True if positive class
    """
    n = len(scores)
    order = sorted(range(n), key=lambda i: scores[i])
    # assign ranks (average for ties)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    n_pos = sum(1 for p in positives if p)
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    sum_ranks_pos = sum(r for r, p in zip(ranks, positives) if p)
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


def summarize_column(name, scores, positives, verdicts_by_row=None):
    auc = roc_auc(scores, positives)
    # simple threshold sweep to find best F1 and best accuracy
    best_f1, best_thr_f1 = -1.0, None
    best_acc, best_thr_acc = -1.0, None
    thrs = sorted(set(scores))
    for thr in thrs:
        tp = fp = fn = tn = 0
        for s, y in zip(scores, positives):
            pred = s >= thr
            if pred and y:
                tp += 1
            elif pred and not y:
                fp += 1
            elif (not pred) and y:
                fn += 1
            else:
                tn += 1
        p = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        r = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0.0
        acc = (tp + tn) / len(scores)
        if f1 > best_f1:
            best_f1, best_thr_f1 = f1, thr
        if acc > best_acc:
            best_acc, best_thr_acc = acc, thr
    return {
        "auc": auc,
        "best_f1": best_f1,
        "best_thr_f1": best_thr_f1,
        "best_acc": best_acc,
        "best_thr_acc": best_thr_acc,
    }


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    json_path, label_path = sys.argv[1], sys.argv[2]
    with open(json_path) as f:
        rows = json.load(f)
    labels = read_labels(label_path)

    # group as chameleon = 'chameleon' OR 'protac'+chameleonic verdict;
    # but for the binary comparison we just use: chameleon | protac = positive
    pos_classes = {"chameleon", "protac"}
    rows_lab = []
    for r in rows:
        lab = labels.get(r["name"])
        if lab is None:
            continue
        rows_lab.append((r, lab))

    print(f"loaded {len(rows_lab)} labelled rows "
          f"(out of {len(rows)} in JSON, "
          f"{len(labels)} in labelled_set)")

    # split out missing
    missing = [r["name"] for r in rows if labels.get(r["name"]) is None]
    if missing:
        print(f"  no label for: {missing}")

    # Gather candidate scalar descriptors
    scalar_keys = []
    for k, v in rows_lab[0][0].items():
        if isinstance(v, (int, float)) and v is not None:
            scalar_keys.append(k)
    # also explicit delta + CI keys; we'll rank them first
    priority = [
        "chameleonic_index", "chameleonic_index_gb",
        "dPSA", "dPSA_gb", "dRg", "dRg_gb", "dIMHB", "dIMHB_gb",
        "psa3d_std", "psa3d_max", "psa3d_min",
        "rg_max", "rg_min", "rg_std",
        "imhb_mean", "imhb_max", "imhb_frac_ge2",
    ]
    ordered = [k for k in priority if k in scalar_keys] + [
        k for k in scalar_keys if k not in priority
    ]

    # build pos flag
    ys = [lab in pos_classes for _, lab in rows_lab]

    print(f"\npositive={sum(ys)} / {len(ys)}")
    print(f"{'descriptor':<32}{'AUC':>8}{'F1*':>8}{'thr(F1)':>10}"
          f"{'acc*':>8}{'thr(acc)':>12}")
    print("-" * 78)
    results = {}
    for k in ordered:
        scores = []
        ys_k = []
        for (r, lab), y in zip(rows_lab, ys):
            v = r.get(k)
            if v is None:
                continue
            scores.append(float(v))
            ys_k.append(y)
        if len(scores) < 4 or sum(ys_k) == 0 or sum(ys_k) == len(ys_k):
            continue
        res = summarize_column(k, scores, ys_k)
        results[k] = res
        print(f"{k:<32}{res['auc']:>8.3f}{res['best_f1']:>8.3f}"
              f"{res['best_thr_f1']:>10.2f}{res['best_acc']:>8.3f}"
              f"{res['best_thr_acc']:>12.2f}")

    # verdict confusion if present
    print("\n--- verdict table (per-molecule) ---")
    print(f"{'name':<22}{'label':<14}{'CI_mmff':>9}{'verdict_mmff':<22}"
          f"{'CI_gb':>8}{'verdict_gb':<22}")
    for r, lab in sorted(rows_lab, key=lambda x: -x[0]["chameleonic_index"]):
        ci = r.get("chameleonic_index")
        v = r.get("chameleonic_verdict", "")
        ci_gb = r.get("chameleonic_index_gb")
        v_gb = r.get("chameleonic_verdict_gb", "")
        ci_s = f"{ci:+.2f}" if ci is not None else ""
        ci_gs = f"{ci_gb:+.2f}" if ci_gb is not None else ""
        print(f"{r['name']:<22}{lab:<14}{ci_s:>9}  {v:<20}"
              f"{ci_gs:>8}  {v_gb:<20}")

    # save
    out = {"descriptors": results, "n_pos": sum(ys), "n": len(ys)}
    out_path = json_path.replace(".json", "_eval.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
