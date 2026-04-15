
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def run():
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    print("START", flush=True)
    results = []
    for i, row in df_user.iterrows():
        name, smiles = row["name"], row["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        tpsa = Descriptors.TPSA(mol)
        mw = Descriptors.MolWt(mol)
        results.append({
            "name": name, "tpsa": tpsa, "mw": mw
        })
    print(pd.DataFrame(results))

if __name__ == "__main__":
    run()
