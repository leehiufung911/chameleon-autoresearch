"""
Cross-reference Chamelogk paper data (SMILES + scores) with existing benchmark.
Build expanded benchmark TSV files.
"""
import csv
import pdfplumber
import re
import os

PROJECT = 'C:/Users/mic23/prototype-embed/chameleon-research-loop'
CHAMELEON = 'C:/Users/mic23/prototype-embed/chameleon'

# ---- 1. Read Chamelogk SMILES from CSV ----
smiles_csv = f'{PROJECT}/papers/chamelogk/jm3c00823_si_002.csv'
csv_compounds = {}
with open(smiles_csv, 'r') as f:
    reader = csv.reader(f, delimiter=';')
    header = next(reader)
    print(f"CSV header: {header}")
    for row in reader:
        if len(row) >= 4 and row[0].strip():
            name = row[0].strip()
            category = row[1].strip()
            subclass = row[2].strip()
            smiles = row[3].strip()
            csv_compounds[name] = {
                'smiles': smiles,
                'category': category,
                'subclass': subclass,
            }
print(f"\nCSV: {len(csv_compounds)} compounds with SMILES")

# ---- 2. Extract Chamelogk scores from PDF Table S2 ----
pdf_path = f'{PROJECT}/papers/chamelogk/jm3c00823_si_001.pdf'
pdf_scores = {}

with pdfplumber.open(pdf_path) as pdf:
    for page in pdf.pages:
        tables = page.extract_tables()
        for table in tables:
            if not table:
                continue
            header_row = table[0]
            # Look for tables with Chamelogk column
            if header_row and 'Chamelogk' in str(header_row):
                chamelogk_idx = None
                name_idx = None
                for i, h in enumerate(header_row):
                    if h and 'Chamelogk' in str(h):
                        chamelogk_idx = i
                    if h and h in ('Compound', 'E3 ligand'):
                        name_idx = i
                if chamelogk_idx is None or name_idx is None:
                    continue
                for row in table[1:]:
                    if not row or not row[name_idx]:
                        continue
                    name = row[name_idx].strip().replace('\n', ' ')
                    try:
                        score = float(row[chamelogk_idx])
                        pdf_scores[name] = score
                    except (ValueError, TypeError):
                        pass

print(f"\nPDF: {len(pdf_scores)} compounds with Chamelogk scores")
print("\nAll Chamelogk scores from PDF:")
for name, score in sorted(pdf_scores.items(), key=lambda x: x[1], reverse=True):
    print(f"  {name:<40} {score:>6.2f}  {'CHAMELEONIC' if score >= 0.6 else 'non-cham'}")

# ---- 3. Match SMILES to scores ----
# Clean up name matching
name_map = {
    'Cyclosporine': 'Cyclosporine',
    'VHL-032 (S,R,S-AHPC HCl)': None,  # separate from PDF Table S3
}

matched = {}
unmatched_scores = []
unmatched_smiles = []

for name, score in pdf_scores.items():
    # Try exact match first
    if name in csv_compounds:
        matched[name] = {**csv_compounds[name], 'chamelogk': score}
    else:
        # Try fuzzy match
        found = False
        for csv_name in csv_compounds:
            if csv_name.lower().replace(' ', '') == name.lower().replace(' ', ''):
                matched[name] = {**csv_compounds[csv_name], 'chamelogk': score}
                found = True
                break
        if not found:
            unmatched_scores.append((name, score))

for name in csv_compounds:
    if name not in matched and not any(name == m for m in matched):
        found = False
        for m_name in matched:
            if name.lower().replace(' ', '') == m_name.lower().replace(' ', ''):
                found = True
                break
        if not found:
            unmatched_smiles.append(name)

print(f"\n\nMatched: {len(matched)} compounds (SMILES + Chamelogk score)")
print(f"Scores without SMILES: {len(unmatched_scores)}")
for name, score in unmatched_scores:
    print(f"  {name}: {score}")
print(f"SMILES without scores: {len(unmatched_smiles)}")
for name in unmatched_smiles:
    print(f"  {name}")

# ---- 4. Read existing benchmark ----
existing = {}
with open(f'{CHAMELEON}/labelled_set.tsv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    for row in reader:
        label, name, smiles = row[0], row[1], row[2]
        existing[name] = {'label': label, 'smiles': smiles}

print(f"\n\nExisting benchmark: {len(existing)} compounds")
labeled = sum(1 for v in existing.values() if v['label'] in ('chameleon', 'nonchameleon'))
print(f"  Labeled: {labeled}")
protacs = sum(1 for v in existing.values() if v['label'] == 'protac')
print(f"  Unlabeled PROTACs: {protacs}")

# ---- 5. Cross-reference ----
print("\n\n=== CROSS-REFERENCE: Existing benchmark vs Chamelogk data ===")
from rdkit import Chem

def canon(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        return Chem.MolToSmiles(mol)
    return smi

# Build canonical SMILES lookup for matched compounds
chamelogk_by_canon = {}
for name, data in matched.items():
    csmi = canon(data['smiles'])
    chamelogk_by_canon[csmi] = (name, data['chamelogk'], data['subclass'])

# Check existing benchmark against Chamelogk
print("\nExisting benchmark molecules with Chamelogk matches:")
overlap_count = 0
for bname, bdata in existing.items():
    bcanon = canon(bdata['smiles'])
    if bcanon in chamelogk_by_canon:
        cname, cscore, csub = chamelogk_by_canon[bcanon]
        clabel = 'chameleonic' if cscore >= 0.6 else 'non-chameleonic'
        agree = "AGREE" if (bdata['label'] == 'chameleon' and cscore >= 0.6) or \
                           (bdata['label'] == 'nonchameleon' and cscore < 0.6) else \
                "PROTAC" if bdata['label'] == 'protac' else "DISAGREE"
        print(f"  {bname:<25} our_label={bdata['label']:<15} Chamelogk={cscore:>5.2f} ({clabel:<16}) {agree}")
        overlap_count += 1

print(f"\nTotal overlaps: {overlap_count}")

# ---- 6. New compounds available to add ----
print("\n\n=== NEW COMPOUNDS AVAILABLE (not in existing benchmark) ===")
existing_canons = {canon(v['smiles']) for v in existing.values()}
new_compounds = []
for name, data in matched.items():
    csmi = canon(data['smiles'])
    if csmi not in existing_canons:
        clabel = 'chameleonic' if data['chamelogk'] >= 0.6 else 'non-chameleonic'
        new_compounds.append({
            'name': name,
            'smiles': data['smiles'],
            'chamelogk': data['chamelogk'],
            'label': clabel,
            'subclass': data['subclass'],
        })
        mol = Chem.MolFromSmiles(data['smiles'])
        mw = sum(a.GetMass() for a in mol.GetAtoms()) if mol else '?'
        print(f"  {name:<40} Chamelogk={data['chamelogk']:>5.2f}  {clabel:<16}  MW~{mw:.0f}  ({data['subclass']})")

print(f"\nTotal new compounds: {len(new_compounds)}")

# ---- 7. Save expanded Chamelogk dataset ----
outpath = f'{PROJECT}/papers/chamelogk/chamelogk_all_matched.tsv'
with open(outpath, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['name', 'smiles', 'chamelogk', 'label', 'class', 'subclass'])
    for name, data in sorted(matched.items(), key=lambda x: x[1]['chamelogk'], reverse=True):
        label = 'chameleonic' if data['chamelogk'] >= 0.6 else 'non-chameleonic'
        w.writerow([name, data['smiles'], data['chamelogk'], label, data['category'], data['subclass']])
print(f"\nSaved {len(matched)} compounds to {outpath}")

# ---- 8. Label summary ----
chameleonic = sum(1 for d in matched.values() if d['chamelogk'] >= 0.6)
noncham = sum(1 for d in matched.values() if d['chamelogk'] < 0.6)
print(f"\nLabel distribution: {chameleonic} chameleonic, {noncham} non-chameleonic")

# PROTACs specifically
protac_data = {n: d for n, d in matched.items() if d['subclass'] == 'PROTAC'}
print(f"PROTACs: {len(protac_data)} total")
for name, data in sorted(protac_data.items(), key=lambda x: x[1]['chamelogk'], reverse=True):
    label = 'chameleonic' if data['chamelogk'] >= 0.6 else 'non-chameleonic'
    print(f"  {name:<30} Chamelogk={data['chamelogk']:>5.2f}  {label}")
