#!/usr/bin/env python
# Ma=0 hydrosphere. run() uses the param.py COMPILED defaults (config_hyd.xml is ignored
# by this path; load_config is unusable from py3 in this build). After the rebuild the
# compiled hydro defaults are nm=1000, checkpoint=20 (VTK every 20). From-scratch + reading
# Transfer_Atm_400 are arranged by file management (stashed hyd_restart_300.bin + _1600),
# and the Ma=0->Ma=10 sweep is stopped by an external watcher.
from pyatom import Hydrosphere
print("=== HYD Ma=0 -> iter 1000 (fresh iter-400 wind, VTK every 20) ===", flush=True)
Hydrosphere().run()
print("=== HYD done ===", flush=True)
