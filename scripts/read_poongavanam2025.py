"""Read Poongavanam 2025 paper and SI PDFs."""
import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import pdfplumber

base = 'C:/Users/mic23/prototype-embed/chameleon-research-loop/papers/Poongavanam2025'

for fname in ['linker-determined-folding-and-hydrophobic-interactions-explain-a-major-difference-in-protac-cell-permeability.pdf', 'SI.pdf']:
    path = f'{base}/{fname}'
    print(f"\n{'='*60}")
    print(f"FILE: {fname}")
    print('='*60)
    with pdfplumber.open(path) as pdf:
        for i, page in enumerate(pdf.pages):
            text = page.extract_text()
            if text:
                print(f"\n--- Page {i+1} ---")
                print(text[:3000])
            tables = page.extract_tables()
            for j, table in enumerate(tables):
                print(f"\n+++ Table {j+1} on page {i+1} +++")
                for row in table:
                    print(row)
