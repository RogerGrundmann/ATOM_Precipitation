#!/bin/bash
# 2026-09-30: is batch 1's byte difference ONLY the -ffast-math fold of `a + alpha_entry*(b - a)` with alpha_entry
# now a literal 1.0? old = HEAD (e79baa9) with that ONE line changed (-O0); run nm 20 from scratch, 1 thread, and
# compare with output_vretb1_a (the full batch 1). Identical => the fold is the whole difference.
set -u; cd "$(dirname "$0")/.."; rm -f python/VFOLD_DONE
W=$(mktemp -d "${TMPDIR:-/tmp}/atom_fold_XXXXXX")
git worktree add --detach "$W" HEAD >/dev/null
( cd "$W" && python3 - <<'PY'
p='atmosphere/SaturationAdjustment.h'; s=open(p).read()
old='''                    const double alpha_entry = ColdCloud::enabled() ? 1.0
                        : 1.0 / (1.0 + std::exp(-(T - m.t_00) / fade_K));'''
assert s.count(old)==1; open(p,'w').write(s.replace(old,'                    const double alpha_entry = 1.0;'))
PY
  make -j8 OPT=-O0 atm > b.log 2>&1 ) && cp -n "$W/cli/atm" cli/fold_atm; git worktree remove --force "$W"; git worktree prune
cd python; [ -e ../cli/fold_atm ] || { echo build failed; touch VFOLD_DONE; exit 1; }
mkdir output_vfold || { echo "output_vfold exists"; touch VFOLD_DONE; exit 1; }
sed "s#output_vconv_new/#output_vfold/#" config_vconv_new.xml > config_vfold.xml
env OMP_NUM_THREADS=1 ../cli/fold_atm config_vfold.xml > vfold.log 2>&1; echo "exit $?"
sed -i 's#output_vfold/#OUT/#' output_vfold/RUN_CONFIG.txt
. ./verify_lib.sh; cmp_dirs vfold vretb1_a "FOLD  old + alpha_entry literal == full batch 1"
touch VFOLD_DONE
