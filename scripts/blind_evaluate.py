"""
blind_evaluate.py — Overseer-only evaluation of agent-generated descriptors.

Reads the agent's iter_N_descriptors.json, joins with benchmark labels AND
holdout labels, computes AUC on both sets, and reports a comparison table.
The agent NEVER sees this output — it goes to logs/iter_N_eval.log.

Usage:
    python blind_evaluate.py <descriptors.json> <labelled_set.tsv> <holdout_set.tsv>
"""
import json
import sys
import os


def read_labels(tsv_path):
    """Read label, name, smiles TSV. Returns {name: label}."""
    labels = {}
    with open(tsv_path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            # Format: label\tname\tsmiles
            labels[parts[1]] = parts[0]
    return labels


def roc_auc(scores, positives):
    """Mann-Whitney U AUC (no sklearn dependency)."""
    n = len(scores)
    if n < 2:
        return float("nan")
    order = sorted(range(n), key=lambda i: scores[i])
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


def evaluate_on_set(descriptors, labels, set_name):
    """Evaluate all numeric descriptor columns against a label set."""
    # Join descriptors with labels
    matched = []
    for d in descriptors:
        name = d.get("name")
        label = labels.get(name)
        if label is not None:
            matched.append((d, label == "chameleon"))

    if len(matched) < 4:
        print(f"  {set_name}: Only {len(matched)} molecules matched — skipping")
        return {}

    n_pos = sum(1 for _, y in matched if y)
    n_neg = len(matched) - n_pos
    print(f"  {set_name}: {len(matched)} molecules (pos={n_pos}, neg={n_neg})")

    # Find all numeric columns (skip 'name', 'smiles', etc.)
    numeric_keys = []
    for k, v in matched[0][0].items():
        if isinstance(v, (int, float)) and k not in ("mw",):
            numeric_keys.append(k)

    if not numeric_keys:
        print(f"  {set_name}: No numeric descriptor columns found")
        return {}

    results = {}
    print(f"  {'Descriptor':<30} {'AUC':>8}")
    print(f"  {'-'*38}")
    for k in numeric_keys:
        scores = []
        ys = []
        for d, y in matched:
            v = d.get(k)
            if v is not None and isinstance(v, (int, float)):
                scores.append(float(v))
                ys.append(y)
        if len(scores) < 4 or sum(ys) == 0 or sum(ys) == len(ys):
            continue
        auc = roc_auc(scores, ys)
        results[k] = auc
        print(f"  {k:<30} {auc:>8.3f}")

    return results


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    desc_path = sys.argv[1]
    bench_path = sys.argv[2]
    holdout_path = sys.argv[3]

    with open(desc_path) as f:
        descriptors = json.load(f)

    print(f"Loaded {len(descriptors)} descriptors from {os.path.basename(desc_path)}")
    print()

    bench_labels = read_labels(bench_path)
    holdout_labels = read_labels(holdout_path)

    print("=== BENCHMARK (49 molecules) ===")
    bench_results = evaluate_on_set(descriptors, bench_labels, "benchmark")

    print()
    print("=== HOLDOUT (29 unseen molecules) ===")
    holdout_results = evaluate_on_set(descriptors, holdout_labels, "holdout")

    # Compare
    print()
    print("=== GENERALIZATION CHECK ===")
    print(f"  {'Descriptor':<30} {'Bench':>8} {'Holdout':>8} {'Gap':>8} {'Flag':>10}")
    print(f"  {'-'*64}")
    for k in bench_results:
        b = bench_results[k]
        h = holdout_results.get(k, float("nan"))
        gap = b - h if not (b != b or h != h) else float("nan")
        flag = "OVERFIT?" if gap > 0.1 else ("GOOD" if gap < 0.05 else "")
        print(f"  {k:<30} {b:>8.3f} {h:>8.3f} {gap:>+8.3f} {flag:>10}")


if __name__ == "__main__":
    main()
