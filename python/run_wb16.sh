#!/bin/bash
# 2026-10-03: wb16 = wb15 + ATM_RH_LAND_EAST=1 (user: try the subtropical east/west asymmetry as a knob). 600 from scratch, cli/atm_ea (-O2),
# 8 threads, runs beside wb15 (cli/atm_ss, same source minus this knob) -- one variable against wb15.
# wb14 land 0-15 / 15-35 / 35-65 / 65-90 (iteration ~220): 354 / 1 / 385 / 7 mm/a (NASA 1653 / 643 / 643 / 295).
# PRE-REGISTERED: subtropical east-side land starts near the wet end (RH ~0.7-0.8 instead of ~0.4-0.5); 15-35 land rain rises from ~0 but stays
# far below NASA -- 50-250 mm/a -- because the free troposphere over subtropical land is warm (env theta_es min 351-355 K against the ocean's 345)
# and the cap limits hot land to the ocean's vapour; the gains sit in India / south China / SE US / E Australia / SE Africa-Brazil; west-coast
# deserts (Sahara west, Atacama, Namib, W Australia, California) stay dry; RISK: Oman / Egypt-Red Sea coast / Arabia east rain (known miss).
# 0-15 and 35-65 land and the ocean unchanged against wb15 (within 3 %).
set -u; cd "$(dirname "$0")"; rm -f WB16_DONE
mkdir output_wb16 || { touch WB16_DONE; exit 1; }
[ -e config_wb16.xml ] && { touch WB16_DONE; exit 1; }
sed "s#output_wb7/#output_wb16/#" config_wb7.xml > config_wb16.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ea $(md5sum < ../cli/atm_ea | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 \
    ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1 ATM_RH_SIGMA_SFC=1 ATM_RH_LAND_EAST=1 ../cli/atm_ea config_wb16.xml > wb16.log 2>&1
t=wb16
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "RH-LAND" $t.log | head -3
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_conv mean" $t.log | tail -1
grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
tail -1 output_$t/convergence.csv
touch WB16_DONE
