#!/bin/bash
# 2026-09-30: [fold test for batch 3: ATM_MICRO_NDIM literal 1.0 in coeff_micro += s*(L/u0 - coeff_micro)] -- adapted from run_vfold2.sh.\n# 2026-09-30: batch 2's only byte difference is signed zeros in PresGradForce at iteration 0 (longal slice). Is it
# the -ffast-math fold of s43*x[1] - s13*x[2] with s43/s13 now literal 1/0 at the phi seam? old = HEAD (786c042)
# with ONLY the two seam lines per model made literal, -O0, nm 20 from scratch, 1 thread, compared with
# output_vretb2_a (the full batch 2). Identical => the fold is the whole difference.
set -u; cd "$(dirname "$0")/.."; rm -f python/VFOLD3_DONE
W=$(mktemp -d "${TMPDIR:-/tmp}/atom_fold3_XXXXXX")
git worktree add --detach "$W" HEAD >/dev/null
( cd "$W" && python3 - <<'PY'
p = 'atmosphere/RHS_Atm_Turb.cpp'; s = open(p).read()
old = """    static const double micro_ndim = [](){
                                           return knob::real(knob::ATM_MICRO_NDIM); }();"""
assert s.count(old) == 1; open(p, 'w').write(s.replace(old, "    static const double micro_ndim = 1.0;"))
PY
  make -j8 OPT=-O0 atm > b.log 2>&1 ) && cp -n "$W/cli/atm" cli/fold3_atm; git worktree remove --force "$W"; git worktree prune
cd python; [ -e ../cli/fold3_atm ] || { echo build failed; touch VFOLD3_DONE; exit 1; }
mkdir output_vfold3 || { echo "output_vfold3 exists"; touch VFOLD3_DONE; exit 1; }
sed "s#output_vconv_new/#output_vfold3/#" config_vconv_new.xml > config_vfold3.xml
env OMP_NUM_THREADS=1 ../cli/fold3_atm config_vfold3.xml > vfold3.log 2>&1; echo "exit $?"
sed -i 's#output_vfold3/#OUT/#' output_vfold3/RUN_CONFIG.txt
. ./verify_lib.sh; cmp_dirs vfold3 vretb3_a "FOLD3  batch 2 + literal micro_ndim == full batch 3"
touch VFOLD3_DONE
