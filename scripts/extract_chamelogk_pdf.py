"""Extract Chamelogk scores from SI PDF using PyPDF2 or pdfplumber."""
import sys
import os
os.environ['PYTHONIOENCODING'] = 'utf-8'
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# Try pdfplumber first, then PyPDF2
pdf_path = 'C:/Users/mic23/prototype-embed/chameleon-research-loop/papers/chamelogk/jm3c00823_si_001.pdf'

try:
    import pdfplumber
    print("Using pdfplumber")
    with pdfplumber.open(pdf_path) as pdf:
        for i, page in enumerate(pdf.pages):
            text = page.extract_text()
            if text:
                print(f"\n=== Page {i+1} ===")
                print(text[:2000])
            tables = page.extract_tables()
            for j, table in enumerate(tables):
                print(f"\n--- Table {j+1} on page {i+1} ---")
                for row in table:
                    print(row)
    sys.exit(0)
except ImportError:
    print("pdfplumber not available")

try:
    import PyPDF2
    print("Using PyPDF2")
    with open(pdf_path, 'rb') as f:
        reader = PyPDF2.PdfReader(f)
        for i, page in enumerate(reader.pages):
            text = page.extract_text()
            if text:
                print(f"\n=== Page {i+1} ===")
                print(text[:2000])
    sys.exit(0)
except ImportError:
    print("PyPDF2 not available")

try:
    import fitz  # pymupdf
    print("Using pymupdf")
    doc = fitz.open(pdf_path)
    for i, page in enumerate(doc):
        text = page.get_text()
        if text:
            print(f"\n=== Page {i+1} ===")
            print(text[:2000])
    sys.exit(0)
except ImportError:
    print("pymupdf not available")

print("\nNo PDF library available. Install one: pip install pdfplumber PyPDF2 pymupdf")
sys.exit(1)
