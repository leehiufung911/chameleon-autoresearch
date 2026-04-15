import os
os.environ["PYTHONUNBUFFERED"] = "1"
import sys
import time
# NO stderr redirect — this is what currently fails in opencode
print("Line 1 - stdout only", flush=True)
time.sleep(2)
print("Line 2 - after 2s", flush=True)
time.sleep(2)
print("Line 3 - after 4s", flush=True)
