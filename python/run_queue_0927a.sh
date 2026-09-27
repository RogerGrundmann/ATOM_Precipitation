#!/bin/bash
# 2026-09-27 queue, phase A: scaling test ALONE on an idle machine, then all four byte checks concurrently
# (17 one-thread arms). Phase B (science arms) is launched separately once the scaling result sets thread counts.
set -u; cd "$(dirname "$0")"; rm -f QUEUE_A_DONE
echo "phase A start $(date +%H:%M)"
bash run_scaling.sh > run_scaling.out 2>&1
echo "scaling done $(date +%H:%M)"
for v in vrks vfx4 vohs vosf; do bash run_verify_$v.sh > run_verify_$v.out 2>&1 & done
wait
echo "byte checks done $(date +%H:%M)"
for v in vrks vfx4 vohs vosf; do echo "== $v"; grep -E "PASS|FAIL|exit [^0]" run_verify_$v.out; done
touch QUEUE_A_DONE
