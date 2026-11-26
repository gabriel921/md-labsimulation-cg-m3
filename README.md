# Aβ42 Coarse-Grained MD Workflow (Martini 3, GROMACS)

End-to-end coarse-grained molecular dynamics workflow using the Martini 3 force field: starts from an atomistic PDB structure, maps it to CG with `martinize2`, runs system preparation, energy minimisation, NVT/NPT equilibration and production MD, then hands off to separate analysis pipelines (coarse-grained analyses, and an optional backmapping stage to all-atom). This script works for Aβ42; to change the system just use a PDB structure from a different protein. `run_setup.sh` counts the CG beads automatically after the mapping, so the chain definition (`ANALYSIS_CHAIN_DEFS=("ChainA:auto")`) only needs editing for multi-chain systems (set `N_CHAINS` accordingly — see [Chains and termini](#chains-and-termini-n_chains-and-martinize2-flags)).

Three scripts drive everything:

- **`run_setup.sh`** — master setup + simulation script (CG mapping → prep → min → NVT → NPT → production).
- **`run_setup_analysis.sh`** — master analysis script for the coarse-grained trajectories, run once production trajectories are ready.
- **`run_setup_backmap.sh`** *(optional)* — backmaps a subset of frames of each replica to all-atom with CG2AT2 and runs all-atom Rg, SASA and DSSP on the result (see [Backmapping](#8-backmapping-to-all-atom-and-all-atom-analyses)).

The simulation stages (NPT + production) can run either as **SLURM jobs** (the original, cluster-based workflow) or **locally**, sequentially, in your own shell, with no SLURM involved. Both modes reuse the exact same `grompp`/`mdrun` commands and the same segment-continuation logic — only the submission/dependency mechanism differs.

Companion repositories: the all-atom (CHARMM36m) and SIRAH 2 versions of this workflow, with the same structure and configuration logic.

---

## 1. Requirements

- GROMACS (tested with 2021.2; both CPU and GPU builds are supported)
- Python 3 with [`vermouth`](https://github.com/marrink-lab/vermouth-martinize) (provides `martinize2`, used for the CG mapping and topology)
- Bash
- The Martini 3 force-field files (see [Force fields and GMXLIB](#6-force-fields-and-gmxlib))
- For SLURM mode: a SLURM cluster with `sbatch`, and either a CPU or GPU partition set up (see [SLURM execution](#5-slurm-execution))
- For the analyses: Python 3 with `numpy`, `scipy`, `matplotlib`, `tqdm`, `mdtraj`, `scikit-learn`, `hdbscan` and `MDAnalysis` (a recent 2.x release, needed to read the GROMACS `.tpr`)
- Optional, for backmapping: [CG2AT2](https://github.com/owenvickery/cg2at) (`pip3 install cg2at2`)
- The analysis scripts have their own details — see the comments at the top of `run_setup_analysis.sh` and the scripts it calls.

---

## 2. Installation

Clone the repository and make sure the expected layout is present:

```
.
├── run_setup.sh
├── run_setup_analysis.sh
├── run_setup_backmap.sh
├── watch_npt_submit_prd.sh
├── ab42_1z0q.pdb
├── mdp/
│   ├── ions_martini3.mdp
│   ├── min_martini3.mdp
│   ├── eq_nvt_martini3.mdp
│   ├── eq_npt_martini3.mdp
│   └── prd_martini3.mdp
└── input/
    ├── simulation/
    │   ├── martinize2.py
    │   ├── make_index_martini3.sh
    │   ├── min_martini3.sh
    │   ├── eq_NVT_martini3.sh
    │   ├── run_segments_local.sh
    │   ├── concat_trj_replica.sh
    │   ├── slurm_MD_eq_NPT_martini3.sh
    │   ├── slurm_MD_eq_NPT_martini3.dat
    │   ├── slurm_MD_prd_martini3.sh
    │   └── slurm_MD_prd_martini3.dat
    ├── analysis/                      (coarse-grained analyses)
    │   ├── prep_nowater.sh
    │   ├── 00_index/index_maker.sh
    │   ├── 01_Rg/        script_Rg_replica.sh, plot_rg.py
    │   ├── 02_SASA/      script_SASA_replica.sh, plot_sasa.py
    │   ├── 04_HDBSCAN/   script_hdbscan_pooled.sh, hdbscan_pooled.py
    │   └── 05_Interactions/  script_intracontacts_pooled.sh, intra_contacts_pooled.py
    ├── analysis_aa/                   (all-atom analyses, used after backmapping)
    │   ├── 00_index/index_maker.sh
    │   ├── 01_Rg/        script_Rg_replica.sh, plot_rg.py
    │   ├── 02_SASA/      script_SASA_replica.sh, plot_sasa.py
    │   └── 03_DSSP/      script_dssp_replica.sh, plot_dssp.py
    └── backmap/
        └── script_backmap_replica.sh
```

Place your input (atomistic) PDB file in the repository root (default filename `ab42_1z0q.pdb`, configurable — see below).

`input/analysis_aa/` holds the analysis scripts of the all-atom repository (same names and layout as its `input/analysis/`). They are only needed for the optional backmapping stage.

No build step is required; the scripts call your GROMACS installation directly.

---

## 3. Configuration

All user-facing configuration lives in the **USER CONFIGURATION** block at the top of `run_setup.sh`. Nothing below that block needs editing for normal use.

### Core settings

| Variable | Meaning |
|---|---|
| `PDB` | Input atomistic PDB filename (must sit next to the script) |
| `SYS` | System name prefix used for output files |
| `N_CHAINS` | Number of chains kept from the input PDB (first model only) |
| `GMX` | Path to the `gmx` executable used for local steps (prep, minimisation, NVT, and — in local mode — NPT/production) |
| `GMXLIB` | Library directory holding the Martini 3 files — see [Force fields and GMXLIB](#6-force-fields-and-gmxlib) |
| `MARTINI_DIR` | Folder with the Martini 3 `.itp` files and `water.gro` |
| `MARTINI_FF_ITP`, `MARTINI_SOLV_ITP`, `MARTINI_ION_ITP` | Martini 3 core, solvent and ion `.itp` file names |
| `SOLVENT_GRO` | Martini 3 water box (inside `MARTINI_DIR`) |
| `MARTINIZE_PY`, `MARTINIZE_FF` | `martinize2` script and force-field name (`martini3001`) |
| `POSRES_FC` | Backbone position-restraint force constant (kJ/mol/nm²), written by `martinize2` |
| `MARTINIZE_EXTRA` | Extra `martinize2` flags (default `-ss C -nt`: all coil, no elastic network — suited to IDPs) |
| `BOX_PADDING`, `BOX_TYPE` | Box parameters |
| `SOLVATE_RADIUS` | `solvate -radius`: vdW radius used for Martini 3 water insertion (nm) |
| `ION_CONC`, `ION_POS`, `ION_NEG`, `ION_SOLVENT_GROUP` | Ionisation parameters (Martini 3 `NA`/`CL`, replacing water `W`) |
| `NTHREADS` | Thread count for local steps |
| `MAXWARN` | `grompp -maxwarn` value |
| `REPLICATES` | Number of independent production replicates (each gets its own `4_prd/repN/`) |

### Chains and termini (`N_CHAINS` and `martinize2` flags)

`N_CHAINS` is the number of chains to keep from the input PDB (the first N); the same `martinize2` options apply to all of them. The termini are controlled through `MARTINIZE_EXTRA`:

```
MARTINIZE_EXTRA="-ss C -nt"
```

By default (no flags) `martinize2` charges both termini (zwitterionic). `-nt` neutralises both termini; `-nter NH2-ter` or `-cter COOH-ter` neutralise only the N- or C-terminus. `-ss C` sets the secondary structure to all coil and, with no `-dssp`, no elastic network is applied — appropriate for intrinsically disordered peptides such as Aβ42. Run `python3 input/simulation/martinize2.py -h` for all options.

> **Warning.** Check the termini choice against your system and the Martini 3 conventions before trusting the defaults.

### Position restraints in NPT (`POSRES_NPT`)

```
POSRES_NPT=yes   # or: no
```

- `yes` — keeps the `define = -DPOSRES` line in the NPT `.mdp`, so the protein backbone is position-restrained (force constant `POSRES_FC`) during NPT equilibration.
- `no` — comments that line out, so NPT runs without restraints.

Only the per-run copy of `eq_npt_martini3.mdp` placed in `3_npt/` is modified; the template in `mdp/` is never touched. With `POSRES_NPT=yes`, the `grompp` call in the NPT scripts needs a reference structure for the restraints (`-r`), e.g. the NVT output `nvt.gro`.

### Execution mode

```
EXECUTION_MODE=slurm   # or: local
```

- `slurm` (default) — NPT and production run as self-requeuing SLURM jobs, with a background watcher chaining them.
- `local` — no SLURM: NPT then each production replicate run sequentially, in the foreground, in the current shell.

See [Local execution](#4-local-execution) and [SLURM execution](#5-slurm-execution) below.

### Step counts and block counts

```
NPT_NSTEPS=50           # nsteps written into the NPT .mdp for this run
NPT_DT=0.02             # dt (ps) written into the NPT .mdp for this run
NPT_SEGMENTS=1          # number of NPT blocks (SLURM Segments= / local loop count)

PROD_NSTEPS=5000        # nsteps written into the production .mdp for this run
PROD_DT=0.02            # dt (ps) written into the production .mdp for this run
PROD_SEGMENTS=3         # number of production blocks per replicate
```

These control:

- The `nsteps` and `dt` values written into the **copies** of `eq_npt_martini3.mdp` / `prd_martini3.mdp` used for the run (the template files under `mdp/` are never modified).
- The number of NPT/production blocks — the SLURM `Segments=` value in `slurm_MD_eq_NPT_martini3.dat` / `slurm_MD_prd_martini3.dat`, or the number of sequential loop iterations in local mode.

One block/segment covers `nsteps × dt` of simulation time and continues directly from the previous block's output structure — this is unchanged in both execution modes. The defaults above are a short end-to-end test with the CG time step of 0.02 ps (20 fs); for real runs raise `nsteps`/segments (e.g. 50000 steps × 0.02 ps = 1 ns per segment).

---

## 4. Local execution

Set `EXECUTION_MODE=local` and run:

```bash
bash run_setup.sh
```

Steps 0–5 (CG mapping and topology, box, solvation, ionisation, minimisation, NVT) are identical to SLURM mode. For steps 6–7, the script calls `input/simulation/run_segments_local.sh` directly — once for NPT, then once per production replicate, one after another, entirely in the foreground. Nothing is submitted to a queue and no background watcher is started.

`run_segments_local.sh` reimplements the same per-segment logic the SLURM scripts use (segment numbering from existing `NNN.gro` files, `tinit` computed as `nsteps × dt × (segment − 1)` read straight from the `.mdp`, `grompp` + `mdrun` continuing from the previous segment) — just without `sbatch`, scratch staging, or self-requeue. It is safe to re-run: any segment whose `.gro` already exists is skipped.

Once production finishes, each replicate's trajectory is concatenated automatically (`concat_trj_replica.sh`), producing `4_prd/repN/traj_all_center_PBC.xtc` — the same file the analysis scripts expect in SLURM mode.

Local mode is intended for workstations or interactive nodes; it uses the `GMX`/`NTHREADS` values from the configuration block for every stage, including NPT and production.

---

## 5. SLURM execution

Set `EXECUTION_MODE=slurm` (the default) and run:

```bash
./run_setup.sh
```

1. NPT equilibration is submitted via `slurm_MD_eq_NPT_martini3.sh`, which self-requeues until `Segments=` (from `slurm_MD_eq_NPT_martini3.dat`) is reached.
2. `watch_npt_submit_prd.sh` runs in the background, waits for NPT's last segment, then submits production for every replicate (`slurm_MD_prd_martini3.sh`, self-requeuing the same way).
3. Once each replicate's production finishes, the watcher concatenates its trajectory automatically.

The step-count/block-count configuration variables described above are applied to the copies of the `.mdp`/`.dat` files placed in `3_npt/` and each `4_prd/repN/` *before* the SLURM scripts are invoked, so the dependency chain, self-requeue logic, and directory layout are completely unaffected.

> **Important:** the analysis (and backmapping) scripts only start after `traj_all_center_PBC.xtc` exists in every replica. The watcher creates it once the **last** production segment of that replica has finished, so wait for all segments to be done (`ls 4_prd/rep*/` and `squeue -u $USER`) before running them.

### Running NPT equilibration on GPU

To run the NPT equilibration on a GPU, edit `input/simulation/slurm_MD_eq_NPT_martini3.dat` and change both the GROMACS version (to a GPU-enabled build) and the partition (to a GPU partition on your cluster). Nothing else in the workflow needs to change.

Monitor with:

```bash
tail -f watch_npt_submit_prd_*.log
squeue -u $USER
```

---

## 6. Force fields and GMXLIB

Martini 3 is not bundled with GROMACS. `GMXLIB` points to an external directory that holds the Martini 3 force-field folder:

```bash
export GMXLIB="/path/to/my/gromacs-topologies"
MARTINI_DIR=$GMXLIB/martini_v300     # your Martini 3 folder name
```

`run_setup.sh` copies the Martini 3 `.itp` files it needs (core force field, solvents, ions) next to the topology in every stage folder (`0_prep/`, `3_npt/`, each `4_prd/repN/`) and makes sure `topol.top` includes them, since `martinize2` only writes the includes needed for the protein.

**Notes:**

- `MARTINI_DIR` must contain the files named in `MARTINI_FF_ITP`, `MARTINI_SOLV_ITP`, `MARTINI_ION_ITP` and `SOLVENT_GRO` (`martini_v3.0.0.itp`, `martini_v3.0.0_solvents_v1.itp`, `martini_v3.0.0_ions_v1.itp`, `water.gro`).
- This repository does **not** ship the Martini 3 force field. You are responsible for obtaining it from the Martini developers ([cgmartini.nl](http://cgmartini.nl)) and complying with their license and citation requirements.
- `run_setup.sh` checks that `GMXLIB`, if set, exists and that the Martini 3 files are present, and stops with a clear error otherwise.

---

## 7. Analysis

Once production trajectories exist (`4_prd/repN/traj_all_center_PBC.xtc`), run:

```bash
./run_setup_analysis.sh
```

> **NOTE:** If you ran the simulation with SLURM, the production script gzips each segment's run-input file, and the analysis scripts look for a plain `.tpr`. Before running the analysis (in particular with `EXECUTION_MODE=ssh`, where the jobs run on other machines), issue `gunzip 001.tpr.gz` inside every replica's folder (`4_prd/rep1`, `4_prd/rep2`, `4_prd/rep3`, ...):
>
> ```bash
> gunzip 4_prd/rep*/001.tpr.gz
> ```
>
> The same applies before `run_setup_backmap.sh` (see [section 8](#8-backmapping-to-all-atom-and-all-atom-analyses)).

See the comments at the top of that script for the individual analyses it runs and the scripts it calls: radius of gyration (`01_Rg`), SASA (`02_SASA`), and the pooled HDBSCAN clustering (`04_HDBSCAN`) and intra-chain contact (`05_Interactions`) analyses. The numbering skips `03` because DSSP is atomistic-only (see the note below).

`run_setup_analysis.sh` has its own `EXECUTION_MODE`, independent of `run_setup.sh`:

- `EXECUTION_MODE=ssh` (default) — dispatches analyses across a pool of machines (`MACHINES=(...)`) over SSH, `MAX_JOBS_PER_MACHINE` at a time per machine. Requires passwordless SSH to each machine.
- `EXECUTION_MODE=local` — no SSH, no machine pool. Everything runs on the current machine instead, `LOCAL_MAX_JOBS` analyses at a time. Set `LOCAL_MAX_JOBS=1` to run them strictly one after another, or raise it to run a few at once — useful if you only have your own workstation.

Each analysis produces a `.png` plot alongside its data — per-replica plots land in `4_prd/repN/analysis/<subfolder>/`, and pooled analyses (HDBSCAN, interactions) land in `4_prd/Analysis/<subfolder>/`.

> **Note on `MACHINES`:** the default hostnames (`bio212`, `bio213`, ...) are specific to this repository's original cluster and won't exist elsewhere. Replace them with your own reachable hostnames, or just use `EXECUTION_MODE=local`.

> **Note on CG analyses:** GROMACS has no default "Backbone" group for CG beads (its own group looks for `N`/`CA`/`C` atoms), so `index_maker.sh` builds it from the Martini 3 backbone bead `BB`. The HDBSCAN clustering uses the `BB` beads. The SASA script gives every bead the same radius (0.264 nm, probe 0.191 nm), so it is a coarse estimate that is useful to compare replicas or conditions, **not** an absolute SASA. DSSP is atomistic-only and is not in `ANALYSES[]`; for secondary structure and absolute SASA use the backmapping stage below.

---

## 8. Backmapping to all-atom and all-atom analyses

`run_setup_backmap.sh` is an optional extra stage. It backmaps a subset of frames of every replica from Martini 3 to all-atom CHARMM36 with [CG2AT2](https://github.com/owenvickery/cg2at) (the fragment-based backmapping tool described in the Martini 3 tutorials) and then runs the **all-atom** Rg, SASA and DSSP analyses on the backmapped trajectory.

It runs in two phases, dispatched across your machines exactly like `run_setup_analysis.sh` (`ssh` or `local`):

1. **Backmap** — one job per replica, in `4_prd/repN/backmap/`: frames are extracted from `traj_all_center_PBC.xtc` (about `NUM_FRAMES` of them, spread uniformly), converted one by one with CG2AT2, and assembled into an all-atom trajectory (`traj_all_center_PBC.xtc`) with its own `001.tpr`.
2. **Analysis** — for every replica that finished phase 1, the all-atom `Rg`, `SASA` and `DSSP` scripts from `input/analysis_aa/` run on that trajectory. The analysis phase starts only when all backmap jobs have finished.

### Requirements

- CG2AT2 installed on every machine used (`pip3 install cg2at2`) — run `cg2at -info` to check that the force-field and fragment names below exist in your installation.
- The all-atom analysis scripts in `input/analysis_aa/`.
- `4_prd/repN/traj_all_center_PBC.xtc` in every replica, and the gunzipped `001.tpr` (see the note in [section 7](#7-analysis)).

### Running

Set the options in the **USER CONFIGURATION** block at the top of `run_setup_backmap.sh` and run it from the system root:

```bash
nohup bash run_setup_backmap.sh > backmap_runner.log 2>&1 &
tail -f backmap_runner.log
```

| Setting | Meaning |
|---|---|
| `EXECUTION_MODE`, `MACHINES`, `MAX_JOBS_PER_MACHINE`, `LOCAL_MAX_JOBS` | Same meaning as in `run_setup_analysis.sh` |
| `CG2AT_BIN` | `cg2at` executable. With `EXECUTION_MODE=ssh` use the **full path** (e.g. `~/.local/bin/cg2at`): remote shells do not load your usual `PATH` |
| `CG2AT_FF`, `CG2AT_FG`, `CG2AT_WATER` | Atomistic force field (`charmm36-cg2at`), Martini 3 fragment library (`martini_3-0_charmm36`) and water model passed to CG2AT2 |
| `CG2AT_MODE` | `de_novo` (minimisation + 5 ps NVT, default) or `none` (minimisation only, faster) |
| `NUM_FRAMES` | Approximate number of frames to backmap per replica |
| `FRAME_JOBS` | CG2AT2 runs in parallel per replica |
| `RUN_PHASES` | `backmap`, `analysis` or `backmap,analysis` (default); can be set as an environment variable |
| `RUN_ANALYSES` | Subset of the all-atom analyses, e.g. `01_Rg,02_SASA` |
| `FORCE_BACKMAP` | `1` redoes the backmapping even if `.backmap_done` exists |

Examples:

```bash
RUN_PHASES=backmap  bash run_setup_backmap.sh     # only backmapping
RUN_PHASES=analysis bash run_setup_backmap.sh     # only the analyses (replicas already backmapped)
```

Start with a small test (e.g. `REPLICATES=(rep1)` and `NUM_FRAMES=3`) before launching everything: converting one frame of the Aβ42 monomer took roughly 10 minutes in our tests.

### Output

```
4_prd/repN/backmap/
├── frames_cg/                  CG frames extracted from the production trajectory
├── cg2at/fNNNNN/               one CG2AT2 run per frame
├── topology/                   topology files used to build the .tpr
├── 001.gro, 001.tpr            all-atom structure + run input (protein only)
├── traj_all_center_PBC.xtc     all-atom trajectory of the backmapped frames
├── index.ndx, system_groups_aa.cfg
├── backmap_<date>.log          backmapping log (one line per frame)
├── .backmap_done               written only when everything succeeded
└── analysis/
    ├── 00_index/               all-atom index
    ├── 01_Rg/  02_SASA/  03_DSSP/    plots and data
```

The dispatch manifest (`backmap_dispatch_<date>.log`) is written in the system root. The workflow is resumable: replicas with `.backmap_done` are skipped, and inside a replica the frames already converted are not redone.

### Notes and limitations

- **No atomistic reference is used on purpose.** CG2AT2 can steer the backmapped structure towards a supplied atomistic structure (`-a`); for an intrinsically disordered peptide this would imprint the starting (e.g. NMR) conformation on every frame and bias Rg, SASA and DSSP. Only the *de novo* modes are supported.
- **Interpret DSSP with care.** De novo backmapped structures have idealised geometry, so the secondary structure is indicative rather than quantitative. Rg and SASA on backmapped frames are more robust (and give an absolute SASA).
- Only **single-chain** systems are handled automatically (the all-atom chain range is generated as `1-N` over the whole protein). For multi-chain systems edit `system_groups_aa.cfg` by hand.
- Backmapping is costly: use a coarse frame selection (`NUM_FRAMES`) rather than backmapping every frame.

---

## Citation

If you use Martini 3, please cite: Souza, P. C. T. et al. Martini 3: a general purpose force field for coarse-grained molecular dynamics. *Nat. Methods* 2021, 18, 382–388. doi: 10.1038/s41592-021-01098-3. Please also follow the citation instructions of `martinize2`/`vermouth` and the Martini website ([cgmartini.nl](http://cgmartini.nl)).

If you use the backmapping stage, please cite CG2AT2: Vickery, O. N.; Stansfeld, P. J. CG2AT2: an Enhanced Fragment-Based Approach for Serial Multi-Scale Molecular Dynamics Simulations. *J. Chem. Theory Comput.* 2021, 17, 6472–6482.
