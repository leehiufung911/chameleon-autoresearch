try:
    import sys, os
    sys.path.append(os.path.abspath("chameleon_local"))
    from chameleon import *
    print("Success")
except Exception as e:
    import traceback
    traceback.print_exc()
