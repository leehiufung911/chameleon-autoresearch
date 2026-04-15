import os
os.environ["PYTHONUNBUFFERED"] = "1"
import sys
# THE FIX: redirect stdout to stderr so Bun's shell doesn't buffer it
sys.stdout = sys.stderr

import time
print("Line 1 - immediate")
time.sleep(2)
print("Line 2 - after 2s sleep")
time.sleep(2)
print("Line 3 - after 4s total")

# Also write to file as backup
with open("experiments/test_stderr_fix_output.txt", "w") as f:
    f.write("File output works too\n")
print("Done - file written")
