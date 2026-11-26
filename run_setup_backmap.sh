#!/bin/bash -e
# =============================================================================
# run_setup_backmap.sh   —  Martini 3 version
#
# Optional extra stage, run AFTER production (and after run_setup_analysis.sh,
# or independently of it): backmaps a subset of frames of every replica from
# Martini 3 to all-atom CHARMM36 with CG2AT2, then runs the ALL-ATOM analyses
# (Rg, SASA, DSSP) on the backmapped trajectory. Useful for absolute SASA and
# for DSSP, which are not meaningful on coarse-grained beads.
#
# Two phases, each dispatched across your machine pool exactly like
# run_setup_analysis.sh (SSH + nohup, or local):
#
#   Phase 1 "backmap"   one job per replica  → 4_prd/repN/backmap/
#                       (input/backmap/script_backmap_replica.sh)
#   Phase 2 "analysis"  one job per replica x analysis, run on the backmapped
#                       all-atom pseudo-replica → 4_prd/repN/backmap/analysis/
#                       (all-atom scripts from input/analysis_aa/)
#
# Choose phases with RUN_PHASES="backmap,analysis" (default), "backmap" or
# "analysis" (e.g. to re-run only the analyses after a backmap finished).
#
# Requirements on every machine used:
#   - CG2AT2 (pip3 install cg2at2  or  conda install -c stansfeld_rg cg2at)
#   - GROMACS (same gmx as the simulation)
#   - the all-atom analysis scripts copied into input/analysis_aa/ (see below)
#
# input/analysis_aa/ must hold the ALL-ATOM repository's analysis scripts, with
# the same layout as its input/analysis/ folder:
#   00_index/index_maker.sh   01_Rg/...   02_SASA/...   03_DSSP/...
# (they are not shipped here because they belong to the all-atom repository).
#
# Run this from the system root (same level as 0_prep/, 4_prd/).
# =============================================================================

set -e

# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# ── EXECUTION MODE (same meaning as in run_setup_analysis.sh) ────────────────
EXECUTION_MODE=ssh                                    # "ssh" or "local"
LOCAL_MAX_JOBS=1                                      # concurrency when EXECUTION_MODE=local
SSH_USER=""                                           # leave empty to use your current user
MACHINES=(bio212 bio213 bio214 bio215)                # pool of machines you're allowed to use THIS run
MAX_JOBS_PER_MACHINE=1                                # backmapping is heavy — keep at 1 unless you know better

GROM_BIN_DEFAULT="/gromacs/gromacs-2021.2/bin/gmx"    # gmx binary — keep = GMX in run_setup.sh

# Per-analysis gmx override (all-atom DSSP is pinned to 2021.x in the all-atom
# repository: do_dssp's argument syntax changed afterwards).
declare -A GROM_BIN_OVERRIDE=(
    ["03_DSSP"]="/gromacs/gromacs-2021.2/bin/gmx"
)

REPLICATES=()                                         # e.g. (rep1 rep2) — empty = auto-discover all repN under 4_prd/
RUN_PHASES="${RUN_PHASES:-backmap,analysis}"          # backmap | analysis | backmap,analysis
DRY_RUN=0                                             # 1 = preview only, dispatch nothing
POLL_INTERVAL=60                                      # seconds between checks for finished jobs

# ── BACKMAPPING SETTINGS (passed to script_backmap_replica.sh) ───────────────
CG2AT_BIN="cg2at"                                     # cg2at executable (path or name on PATH, on the machines used)
CG2AT_FF="charmm36-cg2at"                             # atomistic force field  — see:  cg2at -info
CG2AT_FG="martini_3-0_charmm36"                       # CG fragment library    — see:  cg2at -info
CG2AT_WATER="tip3p"                                   # water model passed to CG2AT2
CG2AT_MODE="de_novo"                                  # de_novo = minimise + 5 ps NVT | none = minimise only (faster)
NUM_FRAMES=50                                         # about how many frames per replica to backmap (spread uniformly)
FRAME_JOBS=1                                          # parallel CG2AT2 runs per replica (per machine)
FORCE_BACKMAP=0                                       # 1 = redo backmapping even if repN/backmap/.backmap_done exists
FORCE_NOWATER_REGEN=0                                 # 1 = regenerate 4_prd/Analysis/nowater.tpr/.gro
FORCE_REINDEX=0                                       # 1 = rebuild the all-atom index_analysis.ndx

# All-atom analyses run on the backmapped trajectory. Same
# "subfolder:script:extra_file1,extra_file2" format as run_setup_analysis.sh;
# master copies live at $SCRIPTS_DIR/analysis_aa/<subfolder>/.
AA_ANALYSES=(
    "01_Rg:script_Rg_replica.sh:plot_rg.py"
    "02_SASA:script_SASA_replica.sh:plot_sasa.py"
    "03_DSSP:script_dssp_replica.sh:plot_dssp.py"
)
RUN_ANALYSES="${RUN_ANALYSES:-}"                      # comma-separated subset of the above, e.g. "01_Rg,03_DSSP"; empty = all

SYSTEM_ROOT="$(pwd)"
SCRIPTS_DIR="$SYSTEM_ROOT/input"
AA_SCRIPTS_DIR="$SCRIPTS_DIR/analysis_aa"
BACKMAP_SCRIPT_DIR="$SCRIPTS_DIR/backmap"
PRD_DIR="$SYSTEM_ROOT/4_prd"
POOLED_DIR="$PRD_DIR/Analysis"                        # nowater.tpr/.gro live here (shared with run_setup_analysis.sh)
# ── END USER CONFIGURATION ───────────────────────────────────────────────────

log() { echo -e "[$(date '+%H:%M:%S')] $*"; }
die() { echo -e "[ERROR] $*" >&2; exit 1; }

case "$EXECUTION_MODE" in
    ssh|local) ;;
    *) die "Unknown EXECUTION_MODE='$EXECUTION_MODE' — must be 'ssh' or 'local'." ;;
esac
if [ "$EXECUTION_MODE" == "local" ]; then
    MACHINES=(LOCAL)
    MAX_JOBS_PER_MACHINE=$LOCAL_MAX_JOBS
fi

[ -d "$PRD_DIR" ] || die "production directory not found: $PRD_DIR (run this from the system root)"
[ "${#MACHINES[@]}" -gt 0 ] || die "MACHINES[] is empty — specify at least one machine."

DO_BACKMAP=0; DO_ANALYSIS=0
IFS=',' read -ra _PH <<< "$RUN_PHASES"
for p in "${_PH[@]}"; do
    case "$p" in
        backmap)  DO_BACKMAP=1 ;;
        analysis) DO_ANALYSIS=1 ;;
        *) die "Unknown phase '$p' in RUN_PHASES — use 'backmap', 'analysis' or 'backmap,analysis'." ;;
    esac
done
[ $((DO_BACKMAP + DO_ANALYSIS)) -gt 0 ] || die "RUN_PHASES is empty."

ssh_host() { if [ -n "$SSH_USER" ]; then echo "${SSH_USER}@$1"; else echo "$1"; fi; }

# ── REMOTE/LOCAL EXECUTION HELPERS (identical to run_setup_analysis.sh) ──────
remote_exec() {
    local host="$1" cmd="$2"
    if [ "$EXECUTION_MODE" == "local" ]; then bash -c "$cmd"; else ssh -o BatchMode=yes "$host" "$cmd"; fi
}
# The job is wrapped in its own redirected subshell in BOTH modes: otherwise the
# backgrounded list keeps ssh's stdout open, ssh does not return until the job
# ends, and the runner blocks (jobs would run one after another).
remote_bg() {
    local host="$1" job_cmd="$2"
    if [ "$EXECUTION_MODE" == "local" ]; then
        bash -c "( $job_cmd ) >/dev/null 2>&1 & echo \$!"
    else
        ssh -o BatchMode=yes "$host" "( $job_cmd ) >/dev/null 2>&1 < /dev/null & echo \$!"
    fi
}
remote_alive() {
    local host="$1" pid="$2"
    if [ "$EXECUTION_MODE" == "local" ]; then kill -0 "$pid" 2>/dev/null; else ssh -o BatchMode=yes "$host" "kill -0 $pid" 2>/dev/null; fi
}

# ── FILTER ANALYSES BY RUN_ANALYSES ───────────────────────────────────────────
if [ -n "$RUN_ANALYSES" ]; then
    IFS=',' read -ra WANTED <<< "$RUN_ANALYSES"
    FILTERED=()
    for entry in "${AA_ANALYSES[@]}"; do
        subfolder="${entry%%:*}"
        for w in "${WANTED[@]}"; do
            if [ "$subfolder" == "$w" ]; then FILTERED+=("$entry"); break; fi
        done
    done
    [ "${#FILTERED[@]}" -gt 0 ] || die "RUN_ANALYSES=\"$RUN_ANALYSES\" matched nothing in AA_ANALYSES[] (known: $(for a in "${AA_ANALYSES[@]}"; do echo -n "${a%%:*} "; done))"
    AA_ANALYSES=("${FILTERED[@]}")
fi

# ── DISCOVER REPLICATES ───────────────────────────────────────────────────────
if [ "${#REPLICATES[@]}" -eq 0 ]; then
    mapfile -t REPLICATES < <(find "$PRD_DIR" -maxdepth 1 -mindepth 1 -type d -name 'rep*' -exec basename {} \; | sort)
fi
[ "${#REPLICATES[@]}" -gt 0 ] || die "no replicate folders found under $PRD_DIR"

# ── PREFLIGHT ─────────────────────────────────────────────────────────────────
if [ "$DO_BACKMAP" -eq 1 ]; then
    [ -f "$BACKMAP_SCRIPT_DIR/script_backmap_replica.sh" ] || die "missing $BACKMAP_SCRIPT_DIR/script_backmap_replica.sh"
    [ -f "$SCRIPTS_DIR/analysis/prep_nowater.sh" ]        || die "missing $SCRIPTS_DIR/analysis/prep_nowater.sh (needed to build the protein-only CG reference)"
fi
if [ "$DO_ANALYSIS" -eq 1 ]; then
    [ -f "$AA_SCRIPTS_DIR/00_index/index_maker.sh" ] || die "missing $AA_SCRIPTS_DIR/00_index/index_maker.sh — copy the ALL-ATOM repository's input/analysis/ scripts into $AA_SCRIPTS_DIR/ (see header)."
    for entry in "${AA_ANALYSES[@]}"; do
        _sub="${entry%%:*}"; _rest="${entry#*:}"; _scr="${_rest%%:*}"; _ext="${_rest#*:}"
        [ "$_ext" == "$_rest" ] && _ext=""
        [ -f "$AA_SCRIPTS_DIR/$_sub/$_scr" ] || die "missing $AA_SCRIPTS_DIR/$_sub/$_scr (copy it from the all-atom repository)"
        if [ -n "$_ext" ]; then
            IFS=',' read -ra _EF <<< "$_ext"
            for _e in "${_EF[@]}"; do [ -f "$AA_SCRIPTS_DIR/$_sub/$_e" ] || die "missing $AA_SCRIPTS_DIR/$_sub/$_e"; done
        fi
    done
fi

log "Replicates        : ${REPLICATES[*]}"
log "Phases            : $([ $DO_BACKMAP -eq 1 ] && echo -n 'backmap ')$([ $DO_ANALYSIS -eq 1 ] && echo -n 'analysis')"
[ "$DO_ANALYSIS" -eq 1 ] && log "AA analyses       : $(for a in "${AA_ANALYSES[@]}"; do echo -n "${a%%:*} "; done)"
[ "$DO_BACKMAP" -eq 1 ]  && log "CG2AT2            : -ff $CG2AT_FF -fg $CG2AT_FG -w $CG2AT_WATER -o $CG2AT_MODE, ~$NUM_FRAMES frames/replica"
echo

MANIFEST="backmap_dispatch_$(date +%Y%m%d_%H%M%S).log"
echo "# phase|replica|analysis|machine|dir|pid|logfile" > "$MANIFEST"

# ── NOWATER PREP (blocking) — protein-only CG reference for frame extraction ─
if [ "$DO_BACKMAP" -eq 1 ]; then
    NOWATER_MACHINE="${MACHINES[0]}"
    nowater_host=$(ssh_host "$NOWATER_MACHINE")
    log "=== Preparing nowater.tpr/.gro (blocking, on $NOWATER_MACHINE) ==="
    mkdir -p "$POOLED_DIR"
    cp -f "$SCRIPTS_DIR/analysis/prep_nowater.sh" "$POOLED_DIR/"
    nowater_cmd="cd '$POOLED_DIR' && PRD_DIR='$PRD_DIR' GROM='$GROM_BIN_DEFAULT' FORCE_REGEN='$FORCE_NOWATER_REGEN' bash prep_nowater.sh"
    if [ "$DRY_RUN" -eq 1 ]; then
        log "  [DRY RUN] would run prep_nowater.sh on $NOWATER_MACHINE"
    else
        remote_exec "$nowater_host" "$nowater_cmd" || die "nowater prep failed on $NOWATER_MACHINE — aborting."
        log "  nowater prep OK."
    fi
    echo
fi

# ── GENERIC QUEUE ENGINE ─────────────────────────────────────────────────────
# A phase fills Q_REP[] / Q_ENTRY[] (and sets Q_KIND to "backmap" or
# "analysis"), then calls run_queue. dispatch_job() handles one job per call.
declare -A MACHINE_LOAD
ACTIVE_MACHINE=(); ACTIVE_PID=(); ACTIVE_REP=(); ACTIVE_SUBFOLDER=(); ACTIVE_LOGFILE=()
Q_REP=(); Q_ENTRY=(); Q_KIND=""; QUEUE_IDX=0; NUM_JOBS=0

dispatch_job() {
    local machine="$1"
    local rep="${Q_REP[$QUEUE_IDX]}"
    local entry="${Q_ENTRY[$QUEUE_IDX]}"
    local job_num=$(( QUEUE_IDX + 1 ))
    QUEUE_IDX=$(( QUEUE_IDX + 1 ))

    local host; host=$(ssh_host "$machine")
    local bm_dir="$PRD_DIR/$rep/backmap"
    local subfolder script extras rest ana_dir idx_dir logfile job_cmd pid

    if [ "$Q_KIND" == "backmap" ]; then
        subfolder="backmap"
        log "=== Backmap job $job_num/$NUM_JOBS — $rep → $machine ==="
        mkdir -p "$bm_dir"
        cp -f "$BACKMAP_SCRIPT_DIR/script_backmap_replica.sh" "$bm_dir/"
        if [ "$FORCE_BACKMAP" -ne 1 ] && [ -f "$bm_dir/.backmap_done" ]; then
            log "  $rep already backmapped (.backmap_done present) — skipping (FORCE_BACKMAP=0)."
            echo "backmap|$rep|backmap|$machine|$bm_dir|SKIPPED|-" >> "$MANIFEST"
            return
        fi
        if [ "$DRY_RUN" -eq 1 ]; then
            log "  [DRY RUN] would run script_backmap_replica.sh for $rep on $machine"
            echo "backmap|$rep|backmap|$machine|$bm_dir|DRYRUN|-" >> "$MANIFEST"
            return
        fi
        logfile="$bm_dir/backmap_$(date +%Y%m%d_%H%M%S).log"
        job_cmd="cd '$bm_dir' && GMX='$GROM_BIN_DEFAULT' CG2AT='$CG2AT_BIN' CG2AT_FF='$CG2AT_FF' CG2AT_FG='$CG2AT_FG' CG2AT_WATER='$CG2AT_WATER' CG2AT_MODE='$CG2AT_MODE' NUM_FRAMES='$NUM_FRAMES' FRAME_JOBS='$FRAME_JOBS' FORCE='$FORCE_BACKMAP' nohup bash script_backmap_replica.sh > '$logfile' 2>&1 < /dev/null"
        ana_dir="$bm_dir"
    else
        subfolder="${entry%%:*}"; rest="${entry#*:}"; script="${rest%%:*}"; extras="${rest#*:}"
        [ "$extras" == "$rest" ] && extras=""
        idx_dir="$bm_dir/analysis/00_index"
        ana_dir="$bm_dir/analysis/$subfolder"
        log "=== AA analysis job $job_num/$NUM_JOBS — $rep, $subfolder → $machine ==="
        mkdir -p "$idx_dir" "$ana_dir"
        cp -f "$AA_SCRIPTS_DIR/00_index/index_maker.sh" "$idx_dir/"
        cp -f "$AA_SCRIPTS_DIR/$subfolder/$script" "$ana_dir/"
        if [ -n "$extras" ]; then
            IFS=',' read -ra EXTRA_FILES <<< "$extras"
            for ef in "${EXTRA_FILES[@]}"; do cp -f "$AA_SCRIPTS_DIR/$subfolder/$ef" "$ana_dir/"; done
        fi
        if [ "$DRY_RUN" -eq 1 ]; then
            log "  [DRY RUN] would build the AA index if needed, then nohup $script on $machine"
            echo "analysis|$rep|$subfolder|$machine|$ana_dir|DRYRUN|-" >> "$MANIFEST"
            return
        fi
        # AA index (foreground, skipped if already built). CFG/GRO/NDX_IN are
        # passed explicitly because this folder sits one level deeper than the
        # usual repN/analysis/00_index/.
        if [ "$FORCE_REINDEX" -eq 1 ] || [ ! -f "$idx_dir/index_analysis.ndx" ]; then
            log "  Building all-atom index on $machine (foreground)..."
            if ! remote_exec "$host" "cd '$idx_dir' && GROM='$GROM_BIN_DEFAULT' CFG='$bm_dir/system_groups_aa.cfg' GRO='$bm_dir/001.gro' NDX_IN='$bm_dir/index.ndx' bash index_maker.sh"; then
                log "  ERROR: index_maker.sh failed for $rep — skipping this job."
                echo "analysis|$rep|$subfolder|$machine|$ana_dir|INDEX_FAILED|-" >> "$MANIFEST"
                return
            fi
        else
            log "  Index already present for $rep — skipping (FORCE_REINDEX=0)."
        fi
        local grom_bin="${GROM_BIN_OVERRIDE[$subfolder]:-$GROM_BIN_DEFAULT}"
        [ -n "${GROM_BIN_OVERRIDE[$subfolder]:-}" ] && log "  Using pinned gmx binary for $subfolder: $grom_bin"
        logfile="$ana_dir/${subfolder}_$(date +%Y%m%d_%H%M%S).log"
        # REP_LABEL: without it the scripts label outputs after the folder two levels up
# ("backmap"), losing the replica name. A hyphen (not an underscore) keeps the
# Rg plot's "<group>_<label>" filename parsing intact.
        job_cmd="cd '$ana_dir' && GROM='$grom_bin' REP_LABEL='${rep}-backmap' nohup bash '$script' > '$logfile' 2>&1 < /dev/null"
    fi

    log "  Dispatching on $machine (nohup, background)..."
    pid=$(remote_bg "$host" "$job_cmd")
    if [ -z "$pid" ]; then
        log "  WARNING: could not confirm PID for $rep/$subfolder on $machine — treating as fire-and-forget (check manually)."
        echo "$Q_KIND|$rep|$subfolder|$machine|$ana_dir|UNKNOWN|$logfile" >> "$MANIFEST"
        return
    fi
    log "  → PID $pid on $machine. Log: $logfile"
    echo "$Q_KIND|$rep|$subfolder|$machine|$ana_dir|$pid|$logfile" >> "$MANIFEST"

    ACTIVE_MACHINE+=("$machine"); ACTIVE_PID+=("$pid"); ACTIVE_REP+=("$rep")
    ACTIVE_SUBFOLDER+=("$subfolder"); ACTIVE_LOGFILE+=("$logfile")
    MACHINE_LOAD["$machine"]=$(( MACHINE_LOAD["$machine"] + 1 ))
}

run_queue() {
    NUM_JOBS=${#Q_REP[@]}
    QUEUE_IDX=0
    for m in "${MACHINES[@]}"; do MACHINE_LOAD["$m"]=0; done
    ACTIVE_MACHINE=(); ACTIVE_PID=(); ACTIVE_REP=(); ACTIVE_SUBFOLDER=(); ACTIVE_LOGFILE=()
    local capacity=$(( ${#MACHINES[@]} * MAX_JOBS_PER_MACHINE ))
    log "Jobs to dispatch  : $NUM_JOBS  (capacity $capacity concurrent: ${#MACHINES[@]} machines x $MAX_JOBS_PER_MACHINE)"

    while [ "$QUEUE_IDX" -lt "$NUM_JOBS" ] || [ "${#ACTIVE_PID[@]}" -gt 0 ]; do
        for machine in "${MACHINES[@]}"; do
            while [ "${MACHINE_LOAD[$machine]}" -lt "$MAX_JOBS_PER_MACHINE" ] && [ "$QUEUE_IDX" -lt "$NUM_JOBS" ]; do
                dispatch_job "$machine"
            done
        done
        if [ "$QUEUE_IDX" -ge "$NUM_JOBS" ] && [ "${#ACTIVE_PID[@]}" -eq 0 ]; then break; fi

        log "  ${#ACTIVE_PID[@]} job(s) running, $(( NUM_JOBS - QUEUE_IDX )) still queued — checking again in ${POLL_INTERVAL}s..."
        sleep "$POLL_INTERVAL"

        local idx machine pid host
        for ((idx=${#ACTIVE_PID[@]}-1; idx>=0; idx--)); do
            machine="${ACTIVE_MACHINE[$idx]}"; pid="${ACTIVE_PID[$idx]}"; host=$(ssh_host "$machine")
            if ! remote_alive "$host" "$pid"; then
                log "  ✓ ${ACTIVE_REP[$idx]}/${ACTIVE_SUBFOLDER[$idx]} finished on $machine (PID $pid). Log: ${ACTIVE_LOGFILE[$idx]}"
                MACHINE_LOAD["$machine"]=$(( MACHINE_LOAD["$machine"] - 1 ))
                unset 'ACTIVE_MACHINE[idx]' 'ACTIVE_PID[idx]' 'ACTIVE_REP[idx]' 'ACTIVE_SUBFOLDER[idx]' 'ACTIVE_LOGFILE[idx]'
            fi
        done
        ACTIVE_MACHINE=("${ACTIVE_MACHINE[@]}"); ACTIVE_PID=("${ACTIVE_PID[@]}"); ACTIVE_REP=("${ACTIVE_REP[@]}")
        ACTIVE_SUBFOLDER=("${ACTIVE_SUBFOLDER[@]}"); ACTIVE_LOGFILE=("${ACTIVE_LOGFILE[@]}")
    done
}

# ── PHASE 1 — BACKMAP (one job per replica) ───────────────────────────────────
if [ "$DO_BACKMAP" -eq 1 ]; then
    log "################ PHASE 1 — backmapping ################"
    Q_KIND="backmap"; Q_REP=(); Q_ENTRY=()
    for rep in "${REPLICATES[@]}"; do Q_REP+=("$rep"); Q_ENTRY+=("backmap"); done
    run_queue
    echo
fi

# ── PHASE 2 — ALL-ATOM ANALYSES (replica x analysis) ──────────────────────────
if [ "$DO_ANALYSIS" -eq 1 ]; then
    log "################ PHASE 2 — all-atom analyses ################"
    Q_KIND="analysis"; Q_REP=(); Q_ENTRY=()
    for rep in "${REPLICATES[@]}"; do
        if [ "$DRY_RUN" -ne 1 ] && [ ! -f "$PRD_DIR/$rep/backmap/.backmap_done" ]; then
            log "WARNING: $rep has no backmap/.backmap_done — backmapping missing or failed (check its backmap_*.log). Skipping its analyses."
            continue
        fi
        for entry in "${AA_ANALYSES[@]}"; do Q_REP+=("$rep"); Q_ENTRY+=("$entry"); done
    done
    if [ "${#Q_REP[@]}" -eq 0 ]; then
        log "No replicate is ready for analysis — nothing to do."
    else
        run_queue
    fi
    echo
fi

log "✅ Done. Manifest: $MANIFEST"
log "   Backmapped trajectories : 4_prd/repN/backmap/traj_all_center_PBC.xtc (+ 001.tpr)"
log "   All-atom analysis output: 4_prd/repN/backmap/analysis/<01_Rg|02_SASA|03_DSSP>/"
