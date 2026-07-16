#!/bin/bash
# Fresh Ma=0 atm -> hyd chain, both to iter 400; hyd uses the atm iter-400 surface wind.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATOM_Precipitation/python || exit 1
TF=output_ATOM_Precipitation/0Ma_smooth_Transfer_Atm_400.vwtp

echo "### STEP 1/2: ATM Ma=0 from scratch -> iter 400  $(date +%H:%M) ###"
rm -f "$TF"                       # ensure the gate checks a fresh file
python3 run_atm.py > atm400.log 2>&1
echo "ATM exit=$?  $(date +%H:%M)"

if [ ! -f "$TF" ]; then
  echo "ERROR: $TF missing after ATM run -> NOT starting HYD."
  echo "Transfer snapshots present:"; ls -1 output_ATOM_Precipitation/0Ma_smooth_Transfer_Atm_*.vwtp 2>/dev/null
  exit 2
fi

echo "### Transfer_Atm_400 present. STEP 2/2: HYD Ma=0 from scratch -> iter 400  $(date +%H:%M) ###"
python3 run_hyd.py > hyd400.log 2>&1
echo "HYD exit=$?  $(date +%H:%M)"
echo "### CHAIN COMPLETE  $(date +%H:%M) ###"
