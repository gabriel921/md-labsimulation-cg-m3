#!/bin/bash -e
# =============================================================================
# min_martini3.sh — energy minimisation of the solvated + ionised Martini 3 system
# (restrained protein: the -DPOSRES define lives in min_CG.mdp, hence -r).
#
# Same step the legacy min_martini3.sh performed, but driven by environment
# variables like the all-atom min script, so run_setup.sh can call it.
# Run with the current directory = 1_min/.
#
# Required environment variables:
#   GMX  NTHREADS  MDP  GRO  TOP
# Output: em.gro (+ em.tpr/.edr/.log/.trr)
# =============================================================================
: "${GMX:?GMX not set}"
: "${NTHREADS:?NTHREADS not set}"
: "${MDP:?MDP not set}"
: "${GRO:?GRO not set}"
: "${TOP:?TOP not set}"

"$GMX" grompp -f "$MDP" -c "$GRO" -r "$GRO" -p "$TOP" -o em.tpr
"$GMX" mdrun -nt "$NTHREADS" -deffnm em
