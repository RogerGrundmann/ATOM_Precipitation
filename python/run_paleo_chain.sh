#!/bin/bash
# Paleo chain: for Ma=100 and Ma=200, run ATM (from scratch -> 400) then, ONLY if the
# ATM run finished cleanly with NO NaN/Inf, the follow-up HYD (from scratch -> 400).
# Sequential (each run gets the full core count). All printouts every 100 iters.
#
# HYD gate — the ATM run must satisfy ALL of:
#   (1) exit code 0 (no abort()),
#   (2) printed its "ATM Ma=<Ma> done" completion marker,
#   (3) NO NaN/Inf in its log (atm prints "non-finite ..." on any non-finite cell,
#       and skips writing the transfer on non-finite surface values), and
#   (4) a freshly written iter-400 surface-wind transfer file.
# If the ATM failed on NaNs/Infs, its HYD is skipped.
set +e
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATOM_Precipitation/python || exit 1
OUT=output_ATOM_Precipitation

for Ma in 100 200; do
  echo "############ [$(date +%H:%M:%S)] ATM Ma=${Ma} from scratch -> iter 400 ############"
  TF=$OUT/${Ma}Ma_smooth_Transfer_Atm_400.vwtp
  rm -f "$TF"                                  # fresh-file gate: existence == written this run
  python3 run_atm_ma.py ${Ma} > atm_${Ma}Ma.log 2>&1
  atm_rc=$?

  nan=$(grep -icE "non-finite|NaN/Inf DETECTED|aborting" atm_${Ma}Ma.log)
  done=$(grep -c "ATM Ma=${Ma} done" atm_${Ma}Ma.log)
  echo "[$(date +%H:%M:%S)] ATM Ma=${Ma}: exit=${atm_rc}  done=${done}  nan/inf_hits=${nan}  transfer=$([ -f "$TF" ] && echo ok || echo MISSING)"

  if [ "${atm_rc}" -ne 0 ] || [ "${done}" -eq 0 ] || [ "${nan}" -gt 0 ] || [ ! -f "$TF" ]; then
    echo "[$(date +%H:%M:%S)] >>> ATM Ma=${Ma} did NOT finish clean (NaN/Inf or incomplete) -> SKIPPING HYD Ma=${Ma}"
    continue
  fi

  echo "############ [$(date +%H:%M:%S)] HYD Ma=${Ma} from scratch -> iter 400 ############"
  python3 run_hyd_ma.py ${Ma} > hyd_${Ma}Ma.log 2>&1
  echo "[$(date +%H:%M:%S)] HYD Ma=${Ma} exit=$?"
done

echo "############ [$(date +%H:%M:%S)] PALEO CHAIN COMPLETE ############"
