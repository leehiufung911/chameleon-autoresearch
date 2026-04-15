
import math
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog("rdApp.*")

def embed_conformers(mol, n_conf=50, seed=42):
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 0
    params.pruneRmsThresh = 0.5
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    return mol, list(cids)

def compute_imhb(mol, conf_id, topo_dist, d_max=2.5, ang_min_deg=120.0):
    conf = mol.GetConformer(conf_id)
    coords = np.asarray(conf.GetPositions())
    donor_pairs = []
    acceptors = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in ("N", "O"):
            acceptors.append(a.GetIdx())
            for n in a.GetNeighbors():
                if n.GetSymbol() == "H":
                    donor_pairs.append((a.GetIdx(), n.GetIdx()))
    if not donor_pairs or not acceptors: return 0
    ang_min_rad = math.radians(ang_min_deg)
    count = 0
    used_donors = set()
    for d_idx, h_idx in donor_pairs:
        if h_idx in used_donors: continue
        for a_idx in acceptors:
            if a_idx == d_idx: continue
            if topo_dist[d_idx, a_idx] < 4: continue
            if np.linalg.norm(coords[h_idx] - coords[a_idx]) > d_max: continue
            v1, v2 = coords[d_idx] - coords[h_idx], coords[a_idx] - coords[h_idx]
            nv1, nv2 = np.linalg.norm(v1), np.linalg.norm(v2)
            if nv1 < 1e-6 or nv2 < 1e-6: continue
            if math.acos(max(-1.0, min(1.0, np.dot(v1, v2) / (nv1 * nv2)))) >= ang_min_rad:
                count += 1
                used_donors.add(h_idx)
                break
    return count

def compute_rg(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    coords = np.asarray(conf.GetPositions())
    masses = np.asarray([a.GetMass() for a in mol.GetAtoms()])
    com = np.average(coords, axis=0, weights=masses)
    sq = np.sum((coords - com) ** 2, axis=1)
    return float(np.sqrt(np.sum(masses * sq) / masses.sum()))

def roc_auc(scores, positives):
    n = len(scores)
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

def get_stats(smiles, n_conf=50):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    mol, cids = embed_conformers(mol, n_conf)
    if not cids: return None
    topo = Chem.GetDistanceMatrix(mol)
    
    # Pre-MMFF
    rgs_pre = [compute_rg(mol, cid) for cid in cids]
    imhbs_pre = [compute_imhb(mol, cid, topo) for cid in cids]
    
    # Post-MMFF
    mol_mmff = Chem.Mol(mol)
    AllChem.MMFFOptimizeMoleculeConfs(mol_mmff, numThreads=0)
    rgs_post = [compute_rg(mol_mmff, cid) for cid in cids]
    imhbs_post = [compute_imhb(mol_mmff, cid, topo) for cid in cids]
    
    return {
        "imhb_pre": np.mean(imhbs_pre),
        "imhb_post": np.mean(imhbs_post),
        "rg_pre_max": np.max(rgs_pre),
        "rg_post_min": np.min(rgs_post),
        "rg_collapse": np.max(rgs_pre) / np.min(rgs_post)
    }

# 1. Test User PROTACs
user_protacs = {
    "protac_1": "O=C(C(N1C(C2=CC=CC(NCCOCCOCCOCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=C2C1=O)=O)CC7)NC7=O",
    "protac_2": "O=C(C(N1C(C2=CC=CC(NC(COCCOCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
    "protac_3": "O=C(C(N1C(C2=CC=CC(OCC(NCCCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O"
}

print("Iteration 2: Testing Hypotheses")
print("-" * 60)
print(f"{'Name':<10} {'IMHB_pre':<10} {'IMHB_post':<10} {'Rg_coll':<10}")
for name, smiles in user_protacs.items():
    res = get_stats(smiles, n_conf=50)
    if res:
        print(f"{name:<10} {res['imhb_pre']:10.4f} {res['imhb_post']:10.4f} {res['rg_collapse']:10.4f}")

# 2. Benchmark Evaluation
benchmark_path = "chameleon_local/labelled_set.tsv"
labels = []
scores_imhb = []
scores_coll = []

print("\nEvaluating Benchmark (n=10 conformers per molecule for speed)...")
with open(benchmark_path) as f:
    header = f.readline().strip().split("\t")
    iL = header.index("label")
    iS = header.index("smiles")
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3: continue
        label_str, smiles = parts[iL], parts[iS]
        res = get_stats(smiles, n_conf=10) # VERY LOW for speed
        if res:
            scores_imhb.append(res['imhb_pre'])
            scores_coll.append(res['rg_collapse'])
            labels.append(label_str == 'chameleon')

auc_imhb = roc_auc(scores_imhb, labels)
auc_coll = roc_auc(scores_coll, labels)
print(f"\nBenchmark AUC (Pre-MMFF IMHB): {auc_imhb:.3f}")
print(f"Benchmark AUC (Rg Collapse Ratio): {auc_coll:.3f}")
