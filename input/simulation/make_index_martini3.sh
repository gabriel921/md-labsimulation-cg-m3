#!/bin/bash -e
# =============================================================================
# make_index_martini3.sh
# Builds the index file for the Martini 3 system (called from STEP 3b of
# run_setup_martini3.sh, inside 0_prep/). No make_ndx menus: the groups are
# written straight from the .gro, so the result does not depend on gmx's
# default group numbering.
#
# Inputs (environment variables, set by run_setup_martini3.sh):
#   GRO       solvated + ionised structure (.gro)
#   NDX       output index file name
#   ION_POS   cation residue name (default NA)
#   ION_NEG   anion  residue name (default CL)
#
# Groups written:
#   System     all beads
#   Protein    every bead that is not water/ion  (use for tc-grps / posres)
#   Solvent    W + ions                          (use for tc-grps)
#   Backbone   beads named BB
# =============================================================================
GRO="${GRO:?GRO not set}"
NDX="${NDX:?NDX not set}"
ION_POS="${ION_POS:-NA}"
ION_NEG="${ION_NEG:-CL}"
[ -f "$GRO" ] || { echo "[ERROR] $GRO not found" >&2; exit 1; }

awk -v ip="$ION_POS" -v in_="$ION_NEG" '
    function flush(name, arr, n,    i) {
        printf "[ %s ]\n", name
        for (i = 1; i <= n; i++) printf "%6d%s", arr[i], (i % 15 == 0 || i == n) ? "\n" : ""
    }
    NR > 2 && NF >= 4 && $0 !~ /^ *-?[0-9.]+ +-?[0-9.]+ +-?[0-9.]+ *$/ {
        i++                                     # atom index = line order
        res = substr($0, 6, 5);  gsub(/ /, "", res)
        atm = substr($0, 11, 5); gsub(/ /, "", atm)
        sys[++ns] = i
        if (res == "W" || res == ip || res == in_) sol[++nw] = i
        else prot[++np] = i
        if (atm == "BB") bb[++nb] = i
    }
    END {
        flush("System", sys, ns)
        flush("Protein", prot, np)
        flush("Solvent", sol, nw)
        flush("Backbone", bb, nb)
        printf "Index: %d beads | Protein %d | Solvent %d | Backbone %d\n", ns, np, nw, nb > "/dev/stderr"
    }
' "$GRO" > "$NDX"

[ -s "$NDX" ] || { echo "[ERROR] $NDX is empty" >&2; exit 1; }
