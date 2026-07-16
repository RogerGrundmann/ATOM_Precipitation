#!/usr/bin/env python
# Single paleo HYD time-slice, from scratch, nm=400, all printouts every 100.
# Reads the atm iter-400 surface wind (<Ma>Ma_smooth_Transfer_Atm_400.vwtp).
# Usage: python run_hyd_ma.py <Ma>
import sys
from pyatom import Hydrosphere
Ma = int(sys.argv[1])
h = Hydrosphere()
h.nm = 400
h.checkpoint = 100          # VTK/MinMax printouts every 100
h.panorama_print = 100      # panorama every 100
print("=== HYD Ma=%d from scratch -> iter 400 (printouts every 100) ===" % Ma, flush=True)
h.run_time_slice(Ma)
print("=== HYD Ma=%d done ===" % Ma, flush=True)
