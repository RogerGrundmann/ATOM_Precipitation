#!/usr/bin/env python
# Continue to nm=400 from the iter-200 checkpoint, on the NEW capped topography.
# moist_phys_start_iter=300 (config default), so moist physics activates around
# iter_n~300 -> the 300..400 window is the real test of the slope cap + trimmed
# filter once SaturationAdjustment / ice scheme / MoistConvection are live.
#
# VTK output enabled so ParaView files are available during/after the run:
#   checkpoint=20     -> radial/zonal/longal slice dumps every 20 iters (~10 sets)
#   panorama_print=50 -> full 3D panorama .vts every 50 iters (~4 dumps over the run)
# Files land in output_path/ as Atmosphere_panorama_*.vts and Atmosphere_*_*.vtk.
from pyatom import Atmosphere

a = Atmosphere()
a.restart_from_iter    = 200      # load atm_restart_200.bin (clean dry-spin-up state)
a.time_start           = 0
a.time_end             = 0        # single Ma=0 slice
a.nm                   = 400      # run iter_n 200 -> 400 (covers coastal u seeding + moist activation)
a.checkpoint_save_iter = -1       # don't overwrite checkpoints during this test
a.checkpoint           = 20       # VTK slice dumps every 20 iters
a.panorama_print       = 50       # full panorama .vts every 50 iters (~700 MB each)
a.run()
