
import math
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors, rdFreeSASA, Descriptors

RDLogger.DisableLog("rdApp.*")

def embed_conformers(mol, n_conf=300, seed=42):
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 0
    params.pruneRmsThresh = 0.5
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    return mol, list(cids)

def polar_atom_mask(mol):
    mask = np.zeros(mol.GetNumAtoms(), dtype=bool)
    for a in mol.GetAtoms():
        if a.GetSymbol() in ("N", "O"):
            mask[a.GetIdx()] = True
        elif a.GetSymbol() == "H":
            if any(n.GetSymbol() in ("N", "O") for n in a.GetNeighbors()):
                mask[a.GetIdx()] = True
    return mask

def compute_3d_psa(mol, conf_id, radii, pmask):
    rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id)
    psa = 0.0
    for atom in mol.GetAtoms():
        if pmask[atom.GetIdx()]:
            psa += float(atom.GetProp("SASA"))
    return psa

def compute_rg(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    coords = np.asarray(conf.GetPositions())
    masses = np.asarray([a.GetMass() for a in mol.GetAtoms()])
    com = np.average(coords, axis=0, weights=masses)
    sq = np.sum((coords - com) ** 2, axis=1)
    return float(np.sqrt(np.sum(masses * sq) / masses.sum()))

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

def get_ci(smiles, n_conf=300):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    mw = Descriptors.MolWt(mol)
    mol, cids = embed_conformers(mol, n_conf)
    if not cids: return None
    radii = rdFreeSASA.classifyAtoms(mol)
    pmask = polar_atom_mask(mol)
    topo = Chem.GetDistanceMatrix(mol)
    psas = [compute_3d_psa(mol, cid, radii, pmask) for cid in cids]
    rgs = [compute_rg(mol, cid) for cid in cids]
    imhbs = [compute_imhb(mol, cid, topo) for cid in cids]
    psa_range = max(psas) - min(psas)
    rg_ratio = max(rgs) / max(min(rgs), 1e-6)
    term_psa = psa_range / math.sqrt(mw)
    term_imhb = np.mean(imhbs)
    term_rg = math.log(rg_ratio)
    ci = 2.0 * term_psa + 1.0 * term_imhb + 3.0 * term_rg
    return ci, term_psa, term_imhb, term_rg, max(rgs)

def manual_auc(labels, scores):
    n = len(labels)
    n_pos = sum(labels)
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0: return 0
    combined = sorted(zip(scores, range(n)), key=lambda x: x[0])
    ranks = [0] * n
    for i, (score, idx) in enumerate(combined):
        ranks[idx] = i + 1
    pos_rank_sum = sum(ranks[i] for i, l in enumerate(labels) if l == 1)
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

print("Testing Pre-MMFF CI hypothesis...")

# 1. Test User PROTACs
user_protacs = {
    "protac_1": "O=C(C(N1C(C2=CC=CC(NCCOCCOCCOCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=C2C1=O)=O)CC7)NC7=O",
    "protac_2": "O=C(C(N1C(C2=CC=CC(NC(COCCOCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O",
    "protac_3": "O=C(C(N1C(C2=CC=CC(OCC(NCCCCC(N(CC3)CCN3C(C=C4)=CC=C4C5=NN(C(NC)=O)[C@@H](C)CC6=CC(OC)=C(OC)C=C65)=O)=O)=C2C1=O)=O)CC7)NC7=O"
}

print("\nUser PROTACs (Pre-MMFF):")
print(f"{'Name':<10} {'CI':<6} {'tPSA':<6} {'tIMHB':<6} {'tRg':<6} {'maxRg':<6}")
for name, smiles in user_protacs.items():
    res = get_ci(smiles)
    if res:
        ci, tp, ti, tr, mr = res
        print(f"{name:<10} {ci:6.2f} {tp:6.2f} {ti:6.2f} {tr:6.2f} {mr:6.2f}")

# 2. Benchmark AUC
benchmark_file = "../chameleon/labelled_set.tsv"
labels, scores = [], []
print("\nRunning benchmark (this may take a minute)...")
with open(benchmark_file) as f:
    lines = f.readlines()
    for line in lines:
        parts = line.strip().split("\t")
        if len(parts) < 3: continue
        label_str, name, smiles = parts[0], parts[1], parts[2]
        if label_str == "label": continue
        ci_res = get_ci(smiles, n_conf=100) # Use 100 for benchmark speed
        if ci_res:
            scores.append(ci_res[0])
            labels.append(1 if label_str == "chameleon" else 0)

auc = manual_auc(labels, scores)
print(f"\nBenchmark AUC (Pre-MMFF, nconf=100): {auc:.3f}")
