#!/bin/bash
# 2026-09-30: batch 2's only byte difference is signed zeros in PresGradForce at iteration 0 (longal slice). Is it
# the -ffast-math fold of s43*x[1] - s13*x[2] with s43/s13 now literal 1/0 at the phi seam? old = HEAD (786c042)
# with ONLY the two seam lines per model made literal, -O0, nm 20 from scratch, 1 thread, compared with
# output_vretb2_a (the full batch 2). Identical => the fold is the whole difference.
set -u; cd "$(dirname "$0")/.."; rm -f python/VFOLD2_DONE
W=$(mktemp -d "${TMPDIR:-/tmp}/atom_fold2_XXXXXX")
git worktree add --detach "$W" HEAD >/dev/null
( cd "$W" && python3 - <<'PY'
for p, knob in [('atmosphere/BC_Atm.h', 'ATM_SEAM_PERIODIC'), ('hydrosphere/BC_Hyd.h', 'HYD_SEAM_PERIODIC')]:
    s = open(p).read()
    old = f'''        static const bool seam_periodic = [](){{
            return knob::on(knob::{knob}); }}();
        const double s43 = seam_periodic ? 1.0 : m.c43;
        const double s13 = seam_periodic ? 0.0 : m.c13;'''
    new = f'''        static const bool seam_periodic = [](){{
            return knob::on(knob::{knob}); }}();
        const double s43 = 1.0;
        const double s13 = 0.0;'''
    assert s.count(old) == 1, p; open(p, 'w').write(s.replace(old, new))
PY
  make -j8 OPT=-O0 atm > b.log 2>&1 ) && cp -n "$W/cli/atm" cli/fold2_atm; git worktree remove --force "$W"; git worktree prune
cd python; [ -e ../cli/fold2_atm ] || { echo build failed; touch VFOLD2_DONE; exit 1; }
mkdir output_vfold2 || { echo "output_vfold2 exists"; touch VFOLD2_DONE; exit 1; }
sed "s#output_vconv_new/#output_vfold2/#" config_vconv_new.xml > config_vfold2.xml
env OMP_NUM_THREADS=1 ../cli/fold2_atm config_vfold2.xml > vfold2.log 2>&1; echo "exit $?"
sed -i 's#output_vfold2/#OUT/#' output_vfold2/RUN_CONFIG.txt
. ./verify_lib.sh; cmp_dirs vfold2 vretb2_a "FOLD2  batch 1 + literal seam == full batch 2"
touch VFOLD2_DONE
