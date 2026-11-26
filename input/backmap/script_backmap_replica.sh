#!/bin/bash -e
# =============================================================================
# script_backmap_replica.sh
#
# Backmaps (Martini 3 -> all-atom CHARMM36) a SUBSET of frames of one replica's
# production trajectory with CG2AT2 (the Martini-endorsed fragment-based tool,
# see the "Backmapping systems with CG2AT" Martini 3 tutorial) and assembles
# them into an all-atom "pseudo-replica" that the ordinary all-atom analysis
# scripts (Rg, SASA, DSSP) can run on unchanged.
#
# Run from inside 4_prd/repN/backmap/ (run_setup_backmap.sh does this).
#
# What it produces in this folder:
#   frames_cg/          CG frames extracted from the production trajectory
#   cg2at/fNNNNN/       one CG2AT2 run per frame (FINAL/final_cg2at_de_novo.pdb)
#   topology/           topol_final.top + force-field/.itp files from CG2AT2
#   001.gro, 001.tpr    all-atom structure + run input (protein only)
#   traj_all_center_PBC.xtc   all-atom trajectory of the backmapped frames
#                       (same file name the analysis scripts always look for)
#   index.ndx           default groups of 001.gro (Protein, Backbone, ...)
#   system_groups_aa.cfg  chain/peptide ranges for the all-atom index_maker.sh
#   .backmap_done       marker, written only when everything above succeeded
#
# IMPORTANT — no atomistic reference (-a) is passed to CG2AT2 on purpose. With
# -a, CG2AT2 steers the de novo structure towards the supplied atomistic
# structure; for an intrinsically disordered peptide that would imprint the
# starting (NMR) conformation on every backmapped frame and bias Rg/SASA/DSSP.
# Only the de novo modes (-o none / -o de_novo) are supported here.
#
# Env vars (all optional; defaults match the 4_prd/repN/backmap/ layout):
#   GMX          gmx binary                         (default: gmx)
#   CG2AT        cg2at executable                   (default: cg2at)
#   CG2AT_FF     atomistic force field (-ff)        (default: charmm36-cg2at)
#   CG2AT_FG     CG fragment library (-fg)          (default: martini_3-0_charmm36)
#   CG2AT_WATER  water model (-w)                   (default: tip3p)
#   CG2AT_MODE   none | de_novo                     (default: de_novo)
#                  none    = fragments + minimisation only (fastest)
#                  de_novo = none + short 5 ps NVT relaxation
#   NUM_FRAMES   approx. number of frames to backmap, spread uniformly over the
#                trajectory (default: 50)
#   FRAME_JOBS   CG2AT2 runs in parallel for this replica (default: 1)
#   MAXWARN      grompp -maxwarn for the analysis .tpr (default: 3; the
#                protein-only system is not neutral, which grompp warns about)
#   FORCE        1 = discard previous frames/results and redo everything
#   REP_DIR      replica folder             (default: ..)
#   CG_TRAJ      CG trajectory   (default: $REP_DIR/traj_all_center_PBC.xtc)
#   CG_TPR       protein-only CG .tpr matching CG_TRAJ  (default: ../../Analysis/nowater.tpr)
#   CG_GRO       protein-only CG .gro matching CG_TRAJ  (default: ../../Analysis/nowater.gro)
#   CG_CFG       CG system_groups.cfg (default: ../../../0_prep/system_groups.cfg)
# =============================================================================

set -e

log() { echo -e "[$(date '+%H:%M:%S')] [backmap] $*"; }
die() { echo -e "[ERROR] [backmap] $*" >&2; exit 1; }

BM_DIR="$(pwd)"
abspath() { ( cd "$(dirname "$1")" 2>/dev/null && echo "$(pwd)/$(basename "$1")" ) || echo "$1"; }

GMX="${GMX:-gmx}"
CG2AT="${CG2AT:-cg2at}"
CG2AT_FF="${CG2AT_FF:-charmm36-cg2at}"
CG2AT_FG="${CG2AT_FG:-martini_3-0_charmm36}"
CG2AT_WATER="${CG2AT_WATER:-tip3p}"
CG2AT_MODE="${CG2AT_MODE:-de_novo}"
NUM_FRAMES="${NUM_FRAMES:-50}"
FRAME_JOBS="${FRAME_JOBS:-1}"
MAXWARN="${MAXWARN:-3}"
FORCE="${FORCE:-0}"
REP_DIR="${REP_DIR:-..}"
CG_TRAJ="$(abspath "${CG_TRAJ:-$REP_DIR/traj_all_center_PBC.xtc}")"
CG_TPR="$(abspath "${CG_TPR:-$REP_DIR/../Analysis/nowater.tpr}")"
CG_GRO="$(abspath "${CG_GRO:-$REP_DIR/../Analysis/nowater.gro}")"
CG_CFG="$(abspath "${CG_CFG:-$REP_DIR/../../0_prep/system_groups.cfg}")"

case "$CG2AT_MODE" in
    none) FINAL_PDB_NAME="final_cg2at_de_novo.pdb" ;;
    de_novo) FINAL_PDB_NAME="final_cg2at_de_novo.pdb" ;;
    *) die "CG2AT_MODE='$CG2AT_MODE' — must be 'none' or 'de_novo' (aligned modes need an atomistic reference and are deliberately not supported, see header)." ;;
esac
[[ "$NUM_FRAMES" =~ ^[1-9][0-9]*$ ]] || die "NUM_FRAMES='$NUM_FRAMES' must be a positive integer."
[[ "$FRAME_JOBS" =~ ^[1-9][0-9]*$ ]] || die "FRAME_JOBS='$FRAME_JOBS' must be a positive integer."

command -v "$GMX"   >/dev/null 2>&1 || die "gmx binary not found: $GMX"
command -v "$CG2AT" >/dev/null 2>&1 || die "cg2at not found: $CG2AT (pip3 install cg2at2, or conda install -c stansfeld_rg cg2at)"
for f in "$CG_TRAJ" "$CG_TPR" "$CG_GRO" "$CG_CFG"; do
    [ -f "$f" ] || die "required file not found: $f"
done

log "==============================================="
log "Backmapping — $(basename "$(cd "$REP_DIR" && pwd)")"
log "==============================================="
log "CG trajectory : $CG_TRAJ"
log "CG reference  : $CG_TPR"
log "CG2AT         : $CG2AT  (-ff $CG2AT_FF -fg $CG2AT_FG -w $CG2AT_WATER -o $CG2AT_MODE)"
log "Target frames : ~$NUM_FRAMES   Parallel CG2AT runs: $FRAME_JOBS"
echo

if [ "$FORCE" -eq 1 ]; then
    log "FORCE=1 — removing previous results."
    rm -rf frames_cg cg2at topology traj_aa.pdb first.pdb 001.gro 001.tpr traj_all_center_PBC.xtc index.ndx system_groups_aa.cfg .backmap_done
fi

# ── 1. TRAJECTORY INFO + FRAME SELECTION ─────────────────────────────────────
chk=$("$GMX" check -f "$CG_TRAJ" 2>&1) || true
NATOMS=$(echo "$chk" | awk '/^# Atoms/ {print $3; exit}')
read -r NFRAMES DT_PS < <(echo "$chk" | awk '$1=="Coords" {print $2, $3; exit}')
[ -n "$NATOMS" ] && [ -n "$NFRAMES" ] && [ -n "$DT_PS" ] || die "could not parse atoms/frames/timestep from 'gmx check -f $CG_TRAJ'."

REF_ATOMS=$(sed -n 2p "$CG_GRO" | tr -d '[:space:]')
if [ "$NATOMS" != "$REF_ATOMS" ]; then
    die "trajectory has $NATOMS atoms but the protein-only reference ($CG_GRO) has $REF_ATOMS — this script expects the protein-only traj_all_center_PBC.xtc (same pairing the pooled analyses use)."
fi

STRIDE=$(( (NFRAMES + NUM_FRAMES - 1) / NUM_FRAMES ))
[ "$STRIDE" -ge 1 ] || STRIDE=1
log "Trajectory: $NFRAMES frames, $NATOMS beads, dt=${DT_PS} ps → stride $STRIDE (about $(( (NFRAMES + STRIDE - 1) / STRIDE )) frames)"

mkdir -p frames_cg cg2at
if [ ! -f frames_cg/frame0.gro ]; then
    log "Extracting CG frames..."
    ( cd frames_cg && echo 0 | "$GMX" trjconv -f "$CG_TRAJ" -s "$CG_TPR" -o frame.gro -sep -skip "$STRIDE" > trjconv.log 2>&1 ) \
        || die "trjconv frame extraction failed (see frames_cg/trjconv.log)."
else
    log "CG frames already extracted — reusing (FORCE=1 to redo)."
fi
mapfile -t FRAMES < <(ls frames_cg/frame[0-9]*.gro | sort -V)
[ "${#FRAMES[@]}" -gt 0 ] || die "no frames extracted into frames_cg/."
log "${#FRAMES[@]} CG frame(s) to backmap."
FRAME_DT_PS=$(awk -v s="$STRIDE" -v d="$DT_PS" 'BEGIN { printf "%.6f", s * d }')
T0=$(sed -n '1s/.*t= *\([-0-9.eE+]*\).*/\1/p' "${FRAMES[0]}")
[ -n "$T0" ] || T0=0

# ── 2. CG2AT2 ON EVERY FRAME (resumable: finished frames are skipped) ────────
convert_frame() {
    local gro="$1" b n tag
    b=$(basename "$gro" .gro); n="${b#frame}"; tag=$(printf 'f%05d' "$n")
    [ -f "$BM_DIR/cg2at/$tag/FINAL/$FINAL_PDB_NAME" ] && return 0
    echo "[$(date '+%H:%M:%S')] [backmap]   $tag started"
    if ( cd "$BM_DIR/cg2at" && "$CG2AT" -c "$BM_DIR/frames_cg/$b.gro" \
          -ff "$CG2AT_FF" -fg "$CG2AT_FG" -w "$CG2AT_WATER" -o "$CG2AT_MODE" \
          -loc "$tag" -silent -gmx "$GMX" > "$tag.log" 2>&1 < /dev/null ); then
        echo "[$(date '+%H:%M:%S')] [backmap]   $tag done"
    else
        echo "[$(date '+%H:%M:%S')] [backmap]   $tag FAILED (see cg2at/$tag.log)"
        return 1
    fi
}
export -f convert_frame
export BM_DIR CG2AT CG2AT_FF CG2AT_FG CG2AT_WATER CG2AT_MODE GMX FINAL_PDB_NAME

log "Running CG2AT2 ($FRAME_JOBS at a time)..."
printf '%s\n' "${FRAMES[@]}" | xargs -P "$FRAME_JOBS" -I{} bash -c 'convert_frame "$1"' _ {} || true

MISSING=()
for gro in "${FRAMES[@]}"; do
    b=$(basename "$gro" .gro); n="${b#frame}"; tag=$(printf 'f%05d' "$n")
    [ -f "cg2at/$tag/FINAL/$FINAL_PDB_NAME" ] || MISSING+=("$tag")
done
if [ "${#MISSING[@]}" -gt 0 ]; then
    die "${#MISSING[@]} frame(s) failed in CG2AT2: ${MISSING[*]:0:10}$([ ${#MISSING[@]} -gt 10 ] && echo ' ...') — see cg2at/<frame>.log and the gromacs-outputs files inside each cg2at/<frame>/ folder. Re-run this script: finished frames are skipped."
fi
log "All ${#FRAMES[@]} frame(s) backmapped."

# ── 3. TOPOLOGY (from the first frame) ───────────────────────────────────────
first_b=$(basename "${FRAMES[0]}" .gro); first_tag=$(printf 'f%05d' "${first_b#frame}")
FIRST_FINAL="$BM_DIR/cg2at/$first_tag/FINAL"
rm -rf topology; mkdir -p topology
[ -f "$FIRST_FINAL/topol_final.top" ] || die "topol_final.top not found in $FIRST_FINAL"
cp "$FIRST_FINAL/topol_final.top" topology/
cp "$FIRST_FINAL"/*.itp topology/ 2>/dev/null || true
cp -r "$FIRST_FINAL"/*.ff topology/ 2>/dev/null || true
# CG2AT2 writes its includes relative to a sibling folder, e.g.
#   #include "../FINAL/charmm36-cg2at.ff/forcefield.itp"
# which does not resolve from topology/. Everything the .top needs was copied
# into topology/ above, so drop the "../FINAL/" prefix (the topology is only
# used to build 001.tpr for analysis).
sed -i 's#\.\./FINAL/##g' topology/topol_final.top

# ── 4. ASSEMBLE THE ALL-ATOM TRAJECTORY ──────────────────────────────────────
# One MODEL per frame. The box is taken from the CG frame (CG2AT2 keeps the CG
# box, but writing it explicitly makes the result deterministic). Orthorhombic
# boxes only (the setup scripts use a cubic box).
log "Assembling all-atom trajectory..."
: > traj_aa.pdb
k=0
for gro in "${FRAMES[@]}"; do
    k=$((k + 1))
    b=$(basename "$gro" .gro); tag=$(printf 'f%05d' "${b#frame}")
    t=$(sed -n '1s/.*t= *\([-0-9.eE+]*\).*/\1/p' "$gro"); [ -n "$t" ] || t=0
    read -r bx by bz _ < <(tail -n1 "$gro")
    awk -v t="$t" -v bx="$bx" -v by="$by" -v bz="$bz" -v n="$k" '
        BEGIN {
            printf "TITLE     backmapped frame t= %.5f step= 0\n", t
            printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", bx*10, by*10, bz*10, 90, 90, 90
            printf "MODEL%9d\n", n
        }
        /^(ATOM|HETATM)/ { print }
        END { print "TER"; print "ENDMDL" }' "cg2at/$tag/FINAL/$FINAL_PDB_NAME" >> traj_aa.pdb
done

awk '/^ENDMDL/ { exit } !/^MODEL/ { print }' traj_aa.pdb > first.pdb
"$GMX" editconf -f first.pdb -o 001.gro > editconf.log 2>&1 || die "editconf failed (see editconf.log)."
AA_ATOMS=$(sed -n 2p 001.gro | tr -d '[:space:]')
log "All-atom system: $AA_ATOMS atoms."

cat > aa_dummy.mdp <<'EOF'
; Dummy .mdp — only used to build 001.tpr (masses/charges/topology) for analysis.
integrator      = md
nsteps          = 0
cutoff-scheme   = Verlet
pbc             = xyz
coulombtype     = reaction-field
rcoulomb        = 0.9
vdwtype         = cut-off
rvdw            = 0.9
EOF

( cd topology && "$GMX" grompp -f "$BM_DIR/aa_dummy.mdp" -c "$BM_DIR/001.gro" -p topol_final.top \
      -o "$BM_DIR/001.tpr" -po "$BM_DIR/aa_dummy_out.mdp" -maxwarn "$MAXWARN" > "$BM_DIR/grompp.log" 2>&1 ) \
    || { echo "----- last lines of grompp.log -----" >&2; tail -n 25 "$BM_DIR/grompp.log" >&2; echo "------------------------------------" >&2
         die "grompp failed for 001.tpr (full output in grompp.log; topology files used are in $BM_DIR/topology/)."; }
[ -f 001.tpr ] || die "001.tpr was not produced."

log "Writing traj_all_center_PBC.xtc ($k frames, t0=${T0} ps, spacing ${FRAME_DT_PS} ps)..."
echo 0 | "$GMX" trjconv -f traj_aa.pdb -s 001.tpr -o traj_all_center_PBC.xtc -t0 "$T0" -timestep "$FRAME_DT_PS" > trjconv_aa.log 2>&1 \
    || die "trjconv failed assembling the xtc (see trjconv_aa.log)."

echo q | "$GMX" make_ndx -f 001.gro -o index.ndx > make_ndx.log 2>&1 || die "make_ndx failed (see make_ndx.log)."

# ── 5. system_groups_aa.cfg FOR THE ALL-ATOM index_maker.sh ──────────────────
# Only single-chain systems are handled: the chain range 1-N covers the whole
# (protein-only) all-atom structure.
mapfile -t CGINFO < <( source "$CG_CFG"; echo "${#CHAIN_DEFS[@]}"; echo "${CHAIN_DEFS[0]%%:*}"; echo "${ALL_CHAINS_NAME:-All_Chains}" )
[ "${CGINFO[0]}" = "1" ] || die "$CG_CFG defines ${CGINFO[0]} chains — this script only handles single-chain systems (per-chain atom ranges differ between CG and all-atom and must be set by hand in system_groups_aa.cfg)."
{
    echo "# Auto-generated by script_backmap_replica.sh on $(date)"
    echo "# All-atom ranges of the backmapped (protein-only) system."
    echo
    echo "BACKBONE_GROUP_NAME=\"Backbone\""
    echo
    echo "CHAIN_DEFS=("
    echo "    \"${CGINFO[1]}:1-${AA_ATOMS}\""
    echo ")"
    echo
    echo "PEPTIDE_DEFS=()"
    echo
    echo "BACKBONE_FOR=(${CGINFO[1]})"
    echo
    echo "ALL_CHAINS_NAME=\"${CGINFO[2]}\""
    echo
    echo "CONCAT_HAS_SOLUTE=0"
    echo "CONCAT_SOLUTE_GROUP=\"Protein\""
    echo "CONCAT_CENTER_OUTPUT_GROUP=\"Protein\""
} > system_groups_aa.cfg

# ── 6. SANITY CHECK + DONE MARKER ────────────────────────────────────────────
chk2=$("$GMX" check -f traj_all_center_PBC.xtc 2>&1) || true
OUT_ATOMS=$(echo "$chk2" | awk '/^# Atoms/ {print $3; exit}')
read -r OUT_FRAMES OUT_DT < <(echo "$chk2" | awk '$1=="Coords" {print $2, $3; exit}')
[ "$OUT_ATOMS" = "$AA_ATOMS" ] || die "xtc has ${OUT_ATOMS:-?} atoms but 001.gro has $AA_ATOMS."
[ "$OUT_FRAMES" = "$k" ] || log "WARNING: xtc has ${OUT_FRAMES:-?} frames, expected $k."
log "xtc: $OUT_FRAMES frame(s), $OUT_ATOMS atoms, timestep ${OUT_DT:-?} ps (expected $FRAME_DT_PS ps)."
awk -v a="${OUT_DT:-0}" -v b="$FRAME_DT_PS" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d > 0.01 * b) }' \
    && log "WARNING: xtc timestep differs from the expected frame spacing — time axes in the plots may be off; check with: $GMX check -f traj_all_center_PBC.xtc"

touch .backmap_done
echo
log "✅ Backmapping finished. All-atom pseudo-replica ready in $BM_DIR"
