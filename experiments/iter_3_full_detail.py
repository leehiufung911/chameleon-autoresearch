
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem

sys.path.append("chameleon_local")
import chameleon

def run():
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    print("START", flush=True)
    results = []
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        print(f"Processing {name}...", flush=True)
        s = chameleon.summarize(name, smiles, n_conf=50, verbose=False)
        results.append({
            "name": name,
            "psa_min": s.psa3d_min,
            "psa_max": s.psa3d_max,
            "psa_range": s.psa3d_max - s.psa3d_min,
            "rg_min": s.rg_min,
            "rg_max": s.rg_max,
            "rg_ratio": s.rg_max / s.rg_min,
            "imhb": s.imhb_mean,
            "ci": s.chameleonic_index
        })

    df_res = pd.DataFrame(results)
    print("\n--- FULL DESCRIPTOR DETAIL ---")
    print(df_res)

if __name__ == "__main__":
    run()
