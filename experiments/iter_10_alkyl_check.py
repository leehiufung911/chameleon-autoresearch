
import sys
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_longest_aliphatic_chain(mol):
    # Find all non-ring aliphatic carbons
    # Aliphatic carbon: [C;R0]
    # We want to find the longest contiguous chain of such carbons.
    
    atoms = []
    for a in mol.GetAtoms():
        if a.GetSymbol() == 'C' and not a.IsInRing():
            # Check if it's aliphatic (only single bonds to other carbons or H)
            # Actually, let's just count contiguous non-ring carbons for simplicity
            atoms.append(a.GetIdx())
            
    if not atoms:
        return 0
        
    # Build a subgraph of these atoms
    adj = {idx: [] for idx in atoms}
    for idx in atoms:
        atom = mol.GetAtomWithIdx(idx)
        for n in atom.GetNeighbors():
            n_idx = n.GetIdx()
            if n_idx in adj:
                adj[idx].append(n_idx)
                
    # Find longest path in this subgraph (it's likely a set of trees/lines)
    def get_max_path(start_node, visited):
        visited.add(start_node)
        max_len = 1
        for neighbor in adj[start_node]:
            if neighbor not in visited:
                max_len = max(max_len, 1 + get_max_path(neighbor, visited.copy()))
        return max_len

    longest = 0
    for node in atoms:
        longest = max(longest, get_max_path(node, set()))
    return longest

def run_experiment():
    df_labels = pd.read_csv("chameleon_local/labelled_set.tsv", sep="\t")
    df_user = pd.read_csv("chameleon_local/user_protacs.tsv", sep="\t", names=["name", "smiles"])
    
    results = []
    for i, row in pd.concat([df_labels, df_user]).iterrows():
        name, smiles = row["name"], row["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: continue
        
        longest_c = get_longest_aliphatic_chain(mol)
        label = row.get("label", "user")
        
        results.append({
            "name": name,
            "label": label,
            "longest_c": longest_c
        })
        
    df_res = pd.DataFrame(results)
    df_res = df_res.sort_values("longest_c", ascending=False)
    
    print("Top 20 molecules by longest aliphatic chain:")
    print(df_res.head(20))
    
    print("\nUser PROTACs:")
    print(df_res[df_res["label"] == "user"])
    
    print("\nChameleons with long chains (> 5):")
    print(df_res[(df_res["label"] == "chameleon") & (df_res["longest_c"] > 5)])

if __name__ == "__main__":
    run_experiment()
