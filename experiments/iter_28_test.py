import sys, os

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout = sys.stderr

import sys

sys.stderr.write("TEST: Starting\n")
sys.stderr.flush()

print("TEST: Print to stdout")
sys.stdout.flush()

with open("experiments/iter_28_test_output.txt", "w") as f:
    f.write("TEST OUTPUT\n")

sys.stderr.write("TEST: Done\n")
sys.stderr.flush()
