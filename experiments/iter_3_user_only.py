
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

sys.path.append("chameleon_local")
import chameleon

def run():
    print("START", flush=True)
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    results = []
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        print(f"Processing {name}...", flush=True)
        s = chameleon.summarize(name, smiles, n_conf=10, verbose=False)
        hbd = rdMolDescriptors.CalcNumHBD(Chem.MolFromSmiles(smiles))
        results.append({
            "name": name, "hbd": hbd, "imhb": s.imhb_mean, 
            "sat_ratio": s.imhb_mean / hbd if hbd > 0 else 1.0,
            "ci": s.chameleonic_index
        })
        print(f"  {name}: SatRatio={results[-1]['sat_ratio']:.2f}", flush=True)

    print("\n--- RESULTS ---", flush=True)
    df_res = pd.DataFrame(results)
    print(df_res)

if __name__ == "__main__":
    run()
