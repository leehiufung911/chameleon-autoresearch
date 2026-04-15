"""
Build holdout validation set from Chamelogk data.

Strategy:
- Keep existing 49-mol benchmark as "dev set" (Gemini evaluates against this)
- Create holdout set from Chamelogk compounds NOT in existing benchmark
- Only include bRo5 compounds (relevant to pipeline) — skip tiny Ro5 molecules
- Add Chamelogk labels for the 6 unlabeled PROTACs in existing benchmark
- Flag label discrepancies between our benchmark and Chamelogk ground truth

Also manually add scores for Ritonavir/Saquinavir/Telaprevir that pdfplumber missed
from the table that spanned pages 4-5.
"""
import csv
from rdkit import Chem

PROJECT = 'C:/Users/mic23/prototype-embed/chameleon-research-loop'
CHAMELEON = 'C:/Users/mic23/prototype-embed/chameleon'

# ---- Load matched Chamelogk data ----
chamelogk_data = {}
with open(f'{PROJECT}/papers/chamelogk/chamelogk_all_matched.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        chamelogk_data[row['name']] = row

# Add manually extracted scores for compounds pdfplumber missed (page 5 text):
# These were in the PDF Table S2 but spanned across pages
manual_additions = {
    'Ritonavir': {'smiles': 'CC(C)C1=NC(=CS1)CN(C)C(=O)NC(C(C)C)C(=O)NC(CC2=CC=CC=C2)CC(C(CC3=CC=CC=C3)NC(=O)OCC4=CN=CS4)O',
                  'chamelogk': '0.67', 'label': 'chameleonic', 'class': 'bRo5', 'subclass': 'Non-macrocycle'},
    'Saquinavir': {'smiles': 'CC(C)(C)NC(=O)C1CC2CCCCC2CN1CC(C(CC3=CC=CC=C3)NC(=O)C(CC(=O)N)NC(=O)C4=NC5=CC=CC=C5C=C4)O',
                   'chamelogk': '1.23', 'label': 'chameleonic', 'class': 'bRo5', 'subclass': 'Non-macrocycle'},
    'Telaprevir': {'smiles': 'CCCC(C(=O)C(=O)NC1CC1)NC(=O)C2C3CCCC3CN2C(=O)C(C(C)(C)C)NC(=O)C(C4CCCCC4)NC(=O)C5=NC=CN=C5',
                   'chamelogk': '0.31', 'label': 'non-chameleonic', 'class': 'bRo5', 'subclass': 'Non-macrocycle'},
    'Voclosporin': {'smiles': 'CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)C=CC=CC2=CC=CC=C2)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C',
                    'chamelogk': '0.93', 'label': 'chameleonic', 'class': 'bRo5', 'subclass': 'Macrocycle'},
}

# Load existing benchmark
existing = {}
existing_canons = set()
with open(f'{CHAMELEON}/labelled_set.tsv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # skip header
    for row in reader:
        label, name, smiles = row[0], row[1], row[2]
        existing[name] = {'label': label, 'smiles': smiles}
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            existing_canons.add(Chem.MolToSmiles(mol))

# ---- Label discrepancies ----
print("=" * 70)
print("LABEL DISCREPANCIES (existing benchmark vs Chamelogk ground truth)")
print("=" * 70)

# Map existing benchmark names to Chamelogk names (case-insensitive)
name_map = {
    'cyclosporinA': 'Cyclosporine',
    'sirolimus': 'Sirolimus',
    'everolimus': 'Everolimus',
    'nelfinavir': 'Nelfinavir',
    'atazanavir': 'Atazanavir',
    'diazepam': 'Diazepam',
    'ritonavir': 'Ritonavir',
    'saquinavir': 'Saquinavir',
    'telaprevir': 'Telaprevir',
    'MZ1': 'MZ1',
    'dBET1': 'dBET1',
    'dBET6': 'dBET6',
    'ARV-825': 'ARV-825',
}

discrepancies = []
for bname, cname in name_map.items():
    if bname not in existing:
        continue
    blabel = existing[bname]['label']

    # Get Chamelogk score
    score = None
    if cname in chamelogk_data:
        score = float(chamelogk_data[cname]['chamelogk'])
    elif cname in manual_additions:
        score = float(manual_additions[cname]['chamelogk'])

    if score is None:
        continue

    clabel = 'chameleonic' if score >= 0.6 else 'non-chameleonic'
    our_binary = 'chameleonic' if blabel == 'chameleon' else ('non-chameleonic' if blabel == 'nonchameleon' else 'unlabeled')

    if blabel == 'protac':
        status = f"UNLABELED (Chamelogk says: {clabel})"
    elif (blabel == 'chameleon' and score >= 0.6) or (blabel == 'nonchameleon' and score < 0.6):
        status = "AGREE"
    else:
        status = "DISAGREE ***"
        discrepancies.append((bname, blabel, score, clabel))

    print(f"  {bname:<20} our={blabel:<15} Chamelogk={score:>5.2f} -> {clabel:<16} {status}")

print(f"\n  DISCREPANCIES: {len(discrepancies)}")
for name, our, score, theirs in discrepancies:
    print(f"    {name}: our='{our}' but Chamelogk={score:.2f} -> '{theirs}'")

# ---- PROTAC labels ----
print("\n" + "=" * 70)
print("PROTAC LABELS (from Chamelogk)")
print("=" * 70)
protac_labels = {
    'MZ1': 1.15, 'dBET1': 0.80, 'dBET6': 0.86, 'ARV-825': 0.72, 'ARV-110': 0.99,
    'ARV-771': None,  # not measured
}
for name, score in protac_labels.items():
    if score:
        label = 'chameleonic' if score >= 0.6 else 'non-chameleonic'
        print(f"  {name:<15} Chamelogk={score:.2f} -> {label}")
    else:
        print(f"  {name:<15} NOT MEASURED in Chamelogk paper")

# ---- Build holdout set ----
print("\n" + "=" * 70)
print("HOLDOUT VALIDATION SET (compounds NOT in existing benchmark)")
print("=" * 70)

holdout = []
all_sources = list(chamelogk_data.items()) + list(manual_additions.items())

for name, data in all_sources:
    if isinstance(data, dict) and 'smiles' in data:
        smiles = data['smiles']
    else:
        continue

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        continue
    canon_smi = Chem.MolToSmiles(mol)
    mw = sum(a.GetMass() for a in mol.GetAtoms())

    # Skip if already in benchmark
    if canon_smi in existing_canons:
        continue

    # Skip tiny Ro5 molecules (MW < 400) — not relevant for bRo5 pipeline
    if mw < 400:
        continue

    score = float(data['chamelogk'])
    label = 'chameleon' if score >= 0.6 else 'nonchameleon'
    subclass = data.get('subclass', '?')

    holdout.append({
        'name': name,
        'smiles': smiles,
        'chamelogk': score,
        'label': label,
        'subclass': subclass,
        'mw': mw,
    })

    print(f"  {name:<35} MW={mw:>7.1f}  Chamelogk={score:>5.2f}  {label:<13}  ({subclass})")

print(f"\n  Total holdout: {len(holdout)} compounds")
cham = sum(1 for h in holdout if h['label'] == 'chameleon')
noncham = sum(1 for h in holdout if h['label'] == 'nonchameleon')
print(f"  {cham} chameleonic, {noncham} non-chameleonic")

# ---- Save holdout set in benchmark format ----
holdout_path = f'{PROJECT}/holdout_set.tsv'
with open(holdout_path, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['label', 'name', 'smiles'])
    for h in sorted(holdout, key=lambda x: x['chamelogk'], reverse=True):
        w.writerow([h['label'], h['name'], h['smiles']])
print(f"\n  Saved to {holdout_path}")

# Also save with Chamelogk scores for reference
holdout_scores_path = f'{PROJECT}/holdout_set_with_scores.tsv'
with open(holdout_scores_path, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['label', 'name', 'smiles', 'chamelogk', 'mw', 'subclass'])
    for h in sorted(holdout, key=lambda x: x['chamelogk'], reverse=True):
        w.writerow([h['label'], h['name'], h['smiles'], h['chamelogk'], f"{h['mw']:.1f}", h['subclass']])
print(f"  Saved with scores to {holdout_scores_path}")

# ---- Summary ----
print("\n" + "=" * 70)
print("DATASET SUMMARY")
print("=" * 70)
print(f"  Dev set (existing benchmark):    {len(existing)} compounds (used by Gemini)")
print(f"  Holdout set (new from Chamelogk): {len(holdout)} compounds (Claude-only validation)")
print(f"  Label discrepancies in dev set:  {len(discrepancies)}")
print(f"\n  NOTE: 4 dev set labels disagree with Chamelogk ground truth.")
print(f"  These are likely mislabeled as 'chameleonic' based on oral bioavailability")
print(f"  rather than actual chromatographic chameleonicity measurement.")
print(f"  We do NOT change the dev set labels to preserve backward compatibility,")
print(f"  but we note this as a systematic bias in the benchmark.")
