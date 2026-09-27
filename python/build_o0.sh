#!/bin/bash
# build_o0.sh <rev|WORKTREE> <name>  --  -O0 binaries for BYTE CHECKS (decided 2026-09-27).
# The tree builds -O2 for science runs, but at -O2 with -ffast-math a change to a physics function can reorder its
# floating-point arithmetic and flip last bits / signed zeros (ATM_PRECIP_UPWIND: 4 S_r cells, -0 vs 0), so a strict
# byte check only proves logic identity when BOTH sides are -O0. This builds <rev> (a commit) or WORKTREE (HEAD plus
# the current uncommitted source diff) at OPT=-O0 in a throw-away git worktree, copies cli/atm and cli/hyd to
# cli/<name>_atm and cli/<name>_hyd, and REMOVES the worktree (three leftover builds in /tmp helped fill the disk).
set -eu
cd "$(dirname "$0")/.."
REV=${1:?usage: build_o0.sh <rev|WORKTREE> <name>}; NAME=${2:?usage: build_o0.sh <rev|WORKTREE> <name>}
for f in cli/${NAME}_atm cli/${NAME}_hyd; do [ -e "$f" ] && { echo "$f exists -- choose another name"; exit 1; }; done
W=$(mktemp -d "${TMPDIR:-/tmp}/atom_o0_XXXXXX")
trap 'git worktree remove --force "$W" >/dev/null 2>&1 || rm -rf "$W"; git worktree prune' EXIT
if [ "$REV" = WORKTREE ]; then
    git worktree add --detach "$W" HEAD >/dev/null
    git diff HEAD -- atmosphere hydrosphere lib cli/*.cpp tinyxml2 Makefile python/param.py > "$W/.wt.diff"
    [ -s "$W/.wt.diff" ] && (cd "$W" && git apply .wt.diff)
else
    git worktree add --detach "$W" "$REV" >/dev/null
fi
(cd "$W" && make -j6 OPT=-O0 atm hyd > build_o0.log 2>&1) || { echo "build failed:"; tail -20 "$W/build_o0.log"; exit 1; }
readelf --debug-dump=info "$W/cli/atm" | grep -m1 DW_AT_producer | grep -q -- ' -O[1-3s]' && { echo "not -O0?"; exit 1; }
cp "$W/cli/atm" cli/${NAME}_atm; cp "$W/cli/hyd" cli/${NAME}_hyd
echo "built -O0: cli/${NAME}_atm cli/${NAME}_hyd  from $REV ($(git rev-parse --short HEAD)$( [ "$REV" = WORKTREE ] && echo ' + working tree'))"
