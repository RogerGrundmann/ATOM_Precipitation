"""
Reformat Thermohalin.h to modern C++ style.
Transformations (applied in order):
  1. Tabs → 4 spaces
  2. static_cast<double> for C-style (double) casts
  3. Remove extra spaces inside brackets .x[ i ] → .x[i]
  4. Remove extra spaces inside parentheses ( expr ) → (expr)
  5. K&R brace style: merge lone { onto preceding line (before any trailing comment)
  6. Add #pragma once + OpenMP header guard at top
"""

import re
import sys

SRC = "/home/roger/SynologyDrive/Cloudstation_Notebook/ATOM_Precipitation/hydrosphere/Thermohalin.h"
DST = SRC  # overwrite in-place (backup made below)

with open(SRC, "r") as f:
    text = f.read()

# --- backup ---
with open(SRC + ".bak", "w") as f:
    f.write(text)

lines = text.splitlines()

# ── 1. Tabs → 4 spaces ────────────────────────────────────────────────────────
lines = [l.expandtabs(4) for l in lines]

# ── 2. static_cast<double> conversions ────────────────────────────────────────
# Pattern A: ( double ) ( expr )  →  static_cast<double>(expr)
#   where expr is everything up to the matching closing paren
#   We handle this with a simple regex that replaces `( double ) (` with
#   `static_cast<double>(` and removes the matching `)`:
#   BUT: we must not remove the wrong `)`.  Since the second `(` is always
#   immediately followed by the expression and then `)`, and the content never
#   contains nested parens that would confuse a simple scan, a regex works fine.
#
# Pattern B: ( double ) i   (bare identifier, only in  d_i = ( double ) i; )

def replace_casts(line):
    # Pattern A: ( double ) ( ... )
    # Replace all occurrences on the line (there may be several chained casts)
    # Strategy: iterate until no more matches
    pat_a = re.compile(r'\( double \) \(([^()]*)\)')
    prev = None
    while prev != line:
        prev = line
        line = pat_a.sub(r'static_cast<double>(\1)', line)

    # Pattern B: ( double ) <identifier>  (not followed by `(`)
    pat_b = re.compile(r'\( double \) (?!\()(\w+)')
    line = pat_b.sub(r'static_cast<double>(\1)', line)

    return line

lines = [replace_casts(l) for l in lines]

# ── 3. Remove spaces inside array brackets: .x[ expr ] → .x[expr] ─────────────
def fix_brackets(line):
    # Remove leading/trailing spaces inside [ ]
    # Repeat until stable (handles nested if any)
    pat = re.compile(r'\[\s+([^\[\]]*?)\s+\]')
    prev = None
    while prev != line:
        prev = line
        line = pat.sub(r'[\1]', line)
    return line

lines = [fix_brackets(l) for l in lines]

# ── 4. Remove extra spaces inside parentheses ─────────────────────────────────
# Only for control-flow keywords and expressions — we remove a single leading/
# trailing space inside any ( ... ) that isn't a cast or macro.
def fix_parens(line):
    # Remove space after ( and before ) — only single spaces, iteratively
    pat_open  = re.compile(r'\( ')
    pat_close = re.compile(r' \)')
    prev = None
    while prev != line:
        prev = line
        line = pat_open.sub('(', line)
        line = pat_close.sub(')', line)
    return line

lines = [fix_parens(l) for l in lines]

# ── 5. K&R brace style ────────────────────────────────────────────────────────
# When line[i] ends with a control-flow closing paren (or function sig `)`),
# and line[i+1] is just whitespace + `{`, merge the `{` before any trailing
# comment on line[i].
#
# Handle trailing `//` comment: insert ` {` before it.

def merge_braces(lines):
    result = []
    skip_next = False
    for i, line in enumerate(lines):
        if skip_next:
            skip_next = False
            continue

        stripped = line.rstrip()

        # Check if the *next* line is a lone brace
        if i + 1 < len(lines):
            next_stripped = lines[i + 1].strip()
            if next_stripped == '{':
                # Does current line end with ) or a comment after )?
                # We want to attach the { to the current line.
                # If there's a // comment, insert { before it.
                comment_match = re.search(r'(//.*)', stripped)
                if comment_match:
                    pos = comment_match.start()
                    stripped = stripped[:pos].rstrip() + ' {  ' + stripped[pos:]
                else:
                    stripped = stripped + ' {'
                result.append(stripped)
                skip_next = True
                continue

        result.append(line.rstrip())
    return result

lines = merge_braces(lines)

# ── 6. Add header at top ───────────────────────────────────────────────────────
header = """\
#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

"""

# Strip existing leading comment block and blank lines, keep it but prepend
# the pragma once + omp guard.  Insert after the opening /* ... */ comment.
# Find the end of the opening block comment.
end_of_header_comment = 0
in_block = False
for i, line in enumerate(lines):
    s = line.strip()
    if s.startswith('/*'):
        in_block = True
    if in_block and '*/' in s:
        end_of_header_comment = i + 1
        break

# Insert header directives after the block comment
lines = (lines[:end_of_header_comment]
         + ['']
         + header.splitlines()
         + lines[end_of_header_comment:])

# ── Write result ───────────────────────────────────────────────────────────────
result_text = '\n'.join(lines) + '\n'

with open(DST, "w") as f:
    f.write(result_text)

print(f"Done. Written to {DST}")
print(f"Backup at {SRC}.bak")
