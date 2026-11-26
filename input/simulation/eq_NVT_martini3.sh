#!/bin/bash -e
# =============================================================================
# eq_NVT_martini3.sh — NVT equilibration of the solvent around a restrained
# Martini 3 protein (eq1_CGPROT.mdp; first step generates velocities).
#
# Same step the legacy eq_nvt_martini3.sh performed, but driven by environment
# variables like the all-atom NVT script. Run with the current directory = 2_nvt/
# (run_setup.sh copies the index file there beforehand).
#
# Required environment variables:
#   GMX  NTHREADS  MDP  GRO  TOP  NDX
# Output: nvt.gro and nvt.cpt (+ nvt.tpr/.edr/.log/.xtc)
# =============================================================================
: "${GMX:?GMX not set}"
: "${NTHREADS:?NTHREADS not set}"
: "${MDP:?MDP not set}"
: "${GRO:?GRO not set}"
: "${TOP:?TOP not set}"
: "${NDX:?NDX not set}"

"$GMX" grompp -f "$MDP" -p "$TOP" -po eq1.mdp -n "$NDX" -c "$GRO" -r "$GRO" -o nvt.tpr
"$GMX" mdrun -nt "$NTHREADS" -deffnm nvt
