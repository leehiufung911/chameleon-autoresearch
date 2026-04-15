
import sys
from unittest.mock import MagicMock
m = MagicMock()
sys.modules["rdkit.ML.Cluster"] = m
sys.modules["rdkit.ML.Cluster.Butina"] = m
import os
sys.path.append(os.getcwd())
try:
    import chameleon_local.chameleon
    print("Import successful")
except Exception as e:
    print(f"Import failed: {e}")
