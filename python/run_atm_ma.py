#!/usr/bin/env python
# Single paleo ATM time-slice, from scratch, nm=400, all printouts every 100.
# Usage: python run_atm_ma.py <Ma>
import sys
from pyatom import Atmosphere
Ma = int(sys.argv[1])
a = Atmosphere()
a.nm = 400
a.checkpoint = 100          # VTK/MinMax printouts every 100
a.panorama_print = 100      # panorama every 100
print("=== ATM Ma=%d from scratch -> iter 400 (printouts every 100) ===" % Ma, flush=True)
a.run_time_slice(Ma)
print("=== ATM Ma=%d done ===" % Ma, flush=True)
