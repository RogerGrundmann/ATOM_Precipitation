# verify_lib.sh -- shared helpers for the both-directions byte checks. Source it; do not run it.
#
# *** ARMS RUN IN PARALLEL, AND THAT IS THE DEFAULT BECAUSE IT IS FREE. ***
# A byte check is N one-thread runs in N separate output directories, reading the input files
# read-only. They are independent processes, so running them concurrently cannot change a single
# byte: this model's non-reproducibility is an OpenMP reduction-order race, and OMP_NUM_THREADS is
# 1 in every arm. What it changes is the wall clock -- 4 arms on 24 cores go from ~4T to ~T.
# run_verify.sh was serial from 2026-09-12 to 2026-09-14 for no reason but the order it was
# written in, at a cost of ~35 min per verification. DO NOT RE-SERIALISE THESE.
#
#   arm <tag> <binary> <config> [VAR=val ...]   launch one arm in the background
#   wait_arms                                   wait for all launched arms, report exit codes
#   cmp_dirs <tagA> <tagB> <label>              byte-compare output_<tagA>/ against output_<tagB>/
#   want_differ <tagA> <tagB> <label>           the CONTROL: these MUST differ, or the check is
#                                               vacuous -- see the note at the bottom.
export OMP_NUM_THREADS=1
_ARM_TAGS=""; _ARM_PIDS=""
arm () {
  local tag=$1 bin=$2 cfg=$3; shift 3
  ( env "$@" "$bin" "$cfg" > "$tag.log" 2>&1; echo $? > ".$tag.exit" ) &
  _ARM_TAGS="$_ARM_TAGS $tag"; _ARM_PIDS="$_ARM_PIDS $!"
  echo "  launched $tag  ($(basename "$bin")  ${*:-clean})  pid $!"
}
wait_arms () {
  echo "=== waiting for arms:$_ARM_TAGS  ($(date +%H:%M:%S))"
  wait $_ARM_PIDS
  local bad=0
  for t in $_ARM_TAGS; do
    local e; e=$(cat ".$t.exit" 2>/dev/null || echo "?")
    echo "  $t exit $e"; [ "$e" = 0 ] || bad=1
    grep -qi "nan\|inf detected" "$t.log" && echo "  !! $t: NaN/Inf reported"
  done
  echo "=== all arms done $(date +%H:%M:%S)"
  return $bad
}
_cmp () {  # $1 $2 -> sets _OK/_BAD
  _OK=0; _BAD=0
  for f in "output_$1"/*; do local b="output_$2/$(basename "$f")"
    if cmp -s "$f" "$b"; then _OK=$((_OK+1)); else echo "    DIFFERS: $(basename "$f")"; _BAD=$((_BAD+1)); fi
  done
}
cmp_dirs   () { _cmp "$1" "$2"; echo "$3: identical $_OK / differing $_BAD  -> $([ $_BAD -eq 0 ] && echo PASS || echo FAIL)"; }
want_differ() { _cmp "$1" "$2" >/dev/null 2>&1; _cmp "$1" "$2" >/dev/null
                echo "$3: identical $_OK / differing $_BAD  -> $([ $_BAD -gt 0 ] && echo PASS || echo 'FAIL -- THE KNOB DOES NOT FIRE AT THIS RUN LENGTH, SO THE PASSES ABOVE ARE VACUOUS')"; }
# want_differ exists because byte-identity ALONE passes a no-op. This tree has already reverted
# one "fix" that was a guard which never fired and whose byte check therefore proved nothing
# (UtilsHyd::writeFile's is_final_result, 2026-09-07). If the two branches of the knob under test
# do not differ on the arms being compared, the verification has measured nothing.
