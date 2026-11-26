# Martini 3 notes — technical and theoretical background for this workflow

Companion to the [README](README.md). The README says **how to run** the workflow; this file explains **what the Martini 3 steps do and why**, how to look at the results, and the caveats to keep in mind when interpreting coarse-grained (CG) data of an intrinsically disordered peptide such as Aβ42.

Contents

1. [What Martini 3 is](#1-what-martini-3-is)
2. [Mapping and topology with martinize2](#2-mapping-and-topology-with-martinize2)
3. [Secondary structure is an input, not an output](#3-secondary-structure-is-an-input-not-an-output)
4. [Martini 3 and disordered proteins](#4-martini-3-and-disordered-proteins)
5. [Looking at the system in VMD](#5-looking-at-the-system-in-vmd)
6. [What the CG analyses actually measure](#6-what-the-cg-analyses-actually-measure)
7. [Backmapping to all-atom](#7-backmapping-to-all-atom)
8. [Time in CG simulations](#8-time-in-cg-simulations)
9. [Troubleshooting](#9-troubleshooting)
10. [References](#10-references)

---

## 1. What Martini 3 is

Martini is a coarse-grained force field in which, roughly, **two to four heavy atoms are represented by one bead**. For the solvent, one bead (`W`) stands for about **four water molecules**. Because far fewer particles interact and the energy landscape is smoother, a time step of about **20 fs** can be used (`dt = 0.02` ps in this workflow), which makes microsecond-scale sampling affordable.

Martini 3 beads come in three sizes — **R**egular, **S**mall and **T**iny — and the interaction between any two bead types is taken from a table of interaction levels. The bead names that appear in this workflow's files are:

| Name | Meaning |
|---|---|
| `BB` | Backbone bead, one per residue |
| `SC1`, `SC2`, ... | Side-chain beads (number and types depend on the residue) |
| `W` | Water bead (about 4 water molecules) |
| `NA`, `CL` | Hydrated sodium / chloride ions |

Practical consequences:

- There is **no explicit hydrogen-bond network** and no individual atoms: backbone hydrogen bonds, `φ/ψ` angles and similar atomistic quantities do not exist in a CG trajectory.
- Residues are recognised by GROMACS as `Protein` because `martinize2` writes standard amino-acid residue names. That is why, unlike SIRAH 2, no extra `residuetypes.dat` link is needed.

---

## 2. Mapping and topology with martinize2

`run_setup.sh` (STEP 0) converts the atomistic PDB into a CG structure and topology with `martinize2` (part of `vermouth`):

```
python3 input/simulation/martinize2.py -f <selected.pdb> -ff martini3001 ... $MARTINIZE_EXTRA
```

What it produces and why it matters:

- **`<SYS>_CG.pdb` / `.gro`** — the CG coordinates (one `BB` bead and some `SC*` beads per residue).
- **`topol.top` + `molecule_0.itp`** — the topology: bonds, angles, dihedrals and, if requested, the position-restraint block (`POSRES`) with force constant `POSRES_FC`. `run_setup.sh` then makes sure `topol.top` also includes the Martini 3 core, solvent and ion `.itp` files.
- **Termini** — controlled in `MARTINIZE_EXTRA`: by default both termini are charged; `-nt` neutralises both; `-nter NH2-ter` / `-cter COOH-ter` neutralise one terminus only. Check this choice against the system you want to model.
- **Elastic network** — none is applied (no `-elastic` flag). For folded proteins an elastic network is normally used to keep the tertiary structure; for an IDP it would be wrong, because it would freeze the very flexibility you want to sample.
- **Position restraints** — used only during equilibration (`POSRES_NPT=yes`), to let the solvent relax around the backbone before production. Restraints are meant for equilibration only: check that your `prd_martini3.mdp` does not define `-DPOSRES`.

After solvation and ionisation, `[ molecules ]` in `topol.top` must match the `.gro`; `run_setup.sh` keeps these counts consistent automatically.

---

## 3. Secondary structure is an input, not an output

In Martini, the **backbone bonded parameters depend on the secondary structure assigned at mapping time** (`-ss` in `martinize2`, or `-dssp` to assign it from the structure). The simulation then keeps the backbone in the conformational basin that those parameters describe.

This workflow uses `-ss C` (everything coil), which is appropriate for a disordered peptide. The consequence is important:

> A CG trajectory with `-ss C` **cannot report helix or β-sheet formation in any quantitative way** — secondary structure is not a property the model is free to explore. Do not run DSSP-like analyses on the CG beads and interpret them as secondary-structure populations.

That is why DSSP is not part of the CG analyses (`run_setup_analysis.sh`). If you need secondary-structure information, backmap frames to all-atom (see [section 7](#7-backmapping-to-all-atom)) and treat the result as indicative.

---

## 4. Martini 3 and disordered proteins

Standard Martini 3 parameters are known to make disordered proteins **too compact**, because protein–water interactions are somewhat too weak. A published correction scales up the protein–water interaction strength to obtain more realistic dimensions of IDPs and multidomain proteins (Thomasen et al., see [references](#10-references)); the effect of solute–solvent interaction scaling on monomeric Aβ42 in Martini 3 was also examined by Smardz et al.

Keep in mind:

- This workflow runs the **standard Martini 3 force field**: the scaling correction is not part of the default setup. If you need it, it has to be added to the topology/force-field files you point `MARTINI_DIR` at, and documented in your methods.
- Rg and SASA from the CG runs will therefore tend to be smaller than experimental or all-atom values. Compare **replicas and conditions with each other** rather than reading the numbers as absolute.

---

## 5. Looking at the system in VMD

`traj_all_center_PBC.xtc` contains only the group set in `ANALYSIS_CONCAT_CENTER_OUTPUT_GROUP` (default `Protein`), so it must be loaded on a **protein-only structure with the same number of beads**. The workflow creates one for the pooled analyses:

```bash
vmd 4_prd/Analysis/nowater.gro 4_prd/rep1/traj_all_center_PBC_VMD.xtc
```

`traj_all_center_PBC_VMD.xtc` is the skip-reduced trajectory made by `concat_trj_replica.sh` for visualisation. A few tips:

- VMD does not know the CG bonds. The Martini website provides a VMD helper script (`cg_bonds`) that draws bonds from the topology; check its usage on [cgmartini.nl](http://cgmartini.nl).
- To check the mapping before running anything, load the atomistic input and the CG output together and confirm every residue was mapped and the bead count matches the `Protein CG beads` line printed by `run_setup.sh`:

  ```bash
  vmd -m ab42_1z0q.pdb 0_prep/ab42_CG.pdb
  ```

- Bead-based representations (`VDW` with a bead-sized radius, or `Licorice` plus bonds) are usually clearer than `Lines`.

---

## 6. What the CG analyses actually measure

| Analysis | What it does on CG data | Caveat |
|---|---|---|
| `01_Rg` | `gmx gyrate` on the `ChainA_BB` group (the `BB` beads) | Mass-weighted with the bead masses of the topology; compact-biased (see section 4) |
| `02_SASA` | `gmx sasa` on `Protein` with a generated `vdwradii.dat`: every bead gets radius **0.264 nm**, probe **0.191 nm** | Uniform-radius estimate — compare replicas/conditions, not an absolute SASA |
| `04_HDBSCAN` | PCA of the `BB` bead positions, then HDBSCAN clustering pooled over all replicas; writes a representative (medoid) structure per cluster | Representatives are CG structures; backmap them if you need all-atom models |
| `05_Interactions` | Counts, per frame, residue pairs with any two beads closer than **6 Å**, intra-chain | The cut-off is a starting point (CG beads are larger than atoms); tune it and report it |

Notes on the numbers:

- The SASA radius follows the usual Martini convention: **0.264 nm ≈ 2^(1/6) × σ/2** for a regular bead with σ = 0.47 nm. Small and tiny beads are smaller (about 0.230 and 0.191 nm); the script uses one value for all beads for simplicity.
- `index_maker.sh` creates the backbone group explicitly because GROMACS' own `Backbone` group looks for atoms named `N`, `CA` and `C`, which do not exist in Martini. `BACKBONE_ATOMS` in that script is `BB`.

---

## 7. Backmapping to all-atom

`run_setup_backmap.sh` uses **CG2AT2**, a fragment-based tool: for each CG bead group it fits pre-built atomistic fragments (from the library `martini_3-0_charmm36`), joins them, and relaxes the result with minimisation (and, in `de_novo` mode, a short 5 ps NVT).

Modes and choices:

- **`-o none`** — fragments plus minimisation only (fastest). **`-o de_novo`** — adds a short NVT relaxation (default here).
- **`-a <atomistic.pdb>` is deliberately not used.** With a reference structure CG2AT2 steers the result towards it, which for an IDP would impose the starting (NMR) conformation on every frame and bias Rg, SASA and DSSP.
- An alternative official tool, **Backward**, also exists for Martini 3 (its mapping files were still labelled beta in the Martini tutorial when this workflow was written).

How to read the results:

- **Rg and SASA** of backmapped frames are reasonably robust and give an *absolute* SASA, which the CG script cannot.
- **DSSP** on de novo structures should be treated as **indicative**: the secondary structure of a backmapped frame is partly set by the fragment geometry, not only by the CG dynamics (and the CG dynamics were run with `-ss C`, see section 3).
- Backmapping is expensive (roughly 10 minutes per frame for the Aβ42 monomer in our tests), so use a coarse frame selection (`NUM_FRAMES`) instead of converting everything.
- Only single-chain systems are handled automatically.

---

## 8. Time in CG simulations

The CG time axis is **not real time**: the smoother energy landscape and the coarser solvent make CG dynamics faster than atomistic ones. Reported times are often given as raw CG time, or with an approximate speed-up factor (about 4 in the Martini 2 convention); if you apply a conversion factor, state it explicitly and check what the Martini 3 literature recommends for your observable. All times in the plots of this workflow are raw CG simulation times.

---

## 9. Troubleshooting

Problems we hit while building and testing this workflow:

| Symptom | Cause and fix |
|---|---|
| `[ERROR] concatenated trajectory not found: .../traj_all_center_PBC.xtc` | Production (or the concatenation) has not finished. The analysis scripts need `traj_all_center_PBC.xtc` in every replica; it is created by `concat_trj_replica.sh` after the last segment. Check `ls -a 4_prd/rep*/` for `.concat_done` and `squeue -u $USER` |
| Analysis or backmap cannot find a `.tpr` | In SLURM mode the production script gzips the `.tpr`. Run `gunzip 4_prd/rep*/001.tpr.gz` before `run_setup_analysis.sh` / `run_setup_backmap.sh` |
| Watcher dies at start with a "PRD dat file not found" message | The watcher must use the file names of your force field (`slurm_MD_prd_martini3.*`); a watcher with hard-coded `_allatom` names will not find them |
| `cg2at not found` when using `EXECUTION_MODE=ssh` | Remote shells do not load your `PATH`: set `CG2AT_BIN` to the full path (e.g. `~/.local/bin/cg2at`) |
| Runner "stops" after `Dispatching on ...` and never reports a PID | The SSH connection stayed open until the job ended. Fixed by running the job in a redirected subshell (`remote_bg` in the runners) |
| `grompp` cannot find `.../charmm36-cg2at.ff/forcefield.itp` while building the backmapped `.tpr` | CG2AT2's topology includes the force field through a relative `../FINAL/...` path. `script_backmap_replica.sh` copies the files into `topology/` and removes that prefix |
| Backmap log looks frozen | With `FRAME_JOBS=1` frames are converted one at a time; the log prints one line per frame (`started` / `done`). Check `cg2at/fNNNNN.log` for the current step |
| `missing .../input/backmap/script_backmap_replica.sh` | The per-replica script has to be copied to `input/backmap/` in the system root |

---

## 10. References

- Souza, P. C. T. et al. Martini 3: a general purpose force field for coarse-grained molecular dynamics. *Nat. Methods* 2021, 18, 382–388. doi: 10.1038/s41592-021-01098-3
- Thomasen, F. E.; Pesce, F.; Roesgaard, M. A.; Tesei, G.; Lindorff-Larsen, K. Improving Martini 3 for Disordered and Multidomain Proteins. *J. Chem. Theory Comput.* 2022, 18, 2033–2041.
- Smardz, P.; Anila, M. M.; Rogowski, P.; Li, M. S.; Różycki, B.; Krupa, P. A Practical Guide to All-Atom and Coarse-Grained Molecular Dynamics Simulations Using Amber and Gromacs: A Case Study of Disulfide-Bond Impact on the Intrinsically Disordered Amyloid Beta. *Int. J. Mol. Sci.* 2024, 25(12), 6698.
- Vickery, O. N.; Stansfeld, P. J. CG2AT2: an Enhanced Fragment-Based Approach for Serial Multi-Scale Molecular Dynamics Simulations. *J. Chem. Theory Comput.* 2021, 17, 6472–6482.
- Martini website and tutorials (mapping with `martinize2`, backmapping with CG2AT2 and Backward): [cgmartini.nl](http://cgmartini.nl)
