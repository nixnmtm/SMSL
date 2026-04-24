# STRUCTURE MECHANICS STATISTICAL LEARNING (SMSL)

Elastic network modeling using **Structure Mechanics Statistical Learning**.

- License: BSD

---

## Introduction

SMSL provides a framework for parameterizing multiscale models for protein
structural analysis. In a typical elastic network model (ENM), springs are
defined between coarse-grained sites and normal mode analysis (NMA) is used
to characterize the vibrational behavior of the system. Many ENMs represent
each residue only at the C-alpha position. SMSL extends this idea by
supporting richer multiscale representations for proteins and other molecules
of interest, including nucleic acids, solvent, and ions.

Residues can be represented using:

- alpha-carbon only (CALPHA)
- alpha-carbon plus sidechain (CASIDE)
- amino N, carboxyl O, and chemically specific sidechain atoms (POLAR)

Fluctuation matching relies heavily on CHARMM  and MDAnalysis.

---

## Please Cite

1. Chu JW, Voth GA. *Biophys J.*, 2006 Mar 1; **90(5)** :1572-82.

2. N. Raj, T. Click, H. Yang, and J.-W. Chu, *Computational and Structural Biotechnology Journal*, 2021, **19**, 5309–5320.

3. N. Raj, T. Click, H. Yang, and J.-W. Chu, *Chemical Science*, 2022, **13**, 3688–3696.  
   https://doi.org/10.1039/D1SC06184D

---

## Installation

### 1. Create a conda environment

`requirements.txt` installs the required packages.

```bash
conda create --name <env_name> --file requirements.txt
conda activate <env_name>
```

### 2. Install SMSL

Download the repository and install it:

```bash
pip install <path_to_repo_or_zip>
```

---

## Development

To run all tests:

```bash
tox
```

To combine coverage data across tox environments:

### Windows

```bash
set PYTEST_ADDOPTS=--cov-append
tox
```

### Linux / macOS

```bash
PYTEST_ADDOPTS=--cov-append tox
```

---

## Fluctuation-Matching Workflow

This workflow relies on a correctly configured CHARMM build capable of
handling large coarse-grained fluctuation-matching calculations.

If larger problem sizes are required, use
`charmmsize_increase_before_configuration.py` before configuring and
compiling CHARMM.

### Workflow Summary

The pipeline proceeds through four main stages:

1. Split the all-atom trajectory into overlapping windows
2. Convert each all-atom window into a coarse-grained model 
3. Convert each coarse-grained model into an ENM
4. Run fluctuation matching independently for each window

---

## Example Pipeline Commands

### 1. Split the trajectory into windows

```bash
fluctmatch splittraj \
    -s ./tests/data/system.tpr \
    -f /path/to/system.xtc \
    -b 1 \
    -e 300000 \
    -w 10000 \
    --type MDA \
    --data ./data \
    --n_jobs 8
```

This creates numbered window directories such as:

```text
data/
├── 1/
│   ├── aa.xtc
│   └── splittraj.log
├── 2/
│   ├── aa.xtc
│   └── splittraj.log
├── 3/
│   ├── aa.xtc
│   └── splittraj.log
...
```

### 2. Convert each AA window to a coarse-grained model

```bash
fluctmatch convert \
    -s ./tests/data/system.tpr \
    -f ./data/1/aa.xtc \
    -o ./data/1/ \
    -p cg \
    -m POLAR \
    -m OM \
    -m ADP \
    -c 48 \
    --write
```

### 3. Convert the CG model to an ENM

```bash
fluctmatch convert \
    -s ./data/1/cg.psf \
    -f ./data/1/cg.dcd \
    -m ENM \
    --rmin 0.0 \
    --rmax 8.0 \
    -p ./data/1/fluctmatch \
    -c 48
```

### 4. Run fluctuation matching

#### Fresh run

```bash
fluctmatch run_fm \
    -s ./data/1/fluctmatch.xplor.psf \
    -f ./data/1/cg.dcd \
    -l ./data/1/charmm.log \
    -o ./data/1/ \
    -e /path/to/charmm \
    -t 310 \
    -c 48 \
    --resid \
    -m 0.001 \
    -p fluctmatch \
    -n 300
```

#### Auto-recovery mode

```bash
fluctmatch run_fm \
    -s ./data/1/fluctmatch.xplor.psf \
    -f ./data/1/cg.dcd \
    -l ./data/1/charmm.log \
    -o ./data/1/ \
    -e /path/to/charmm \
    -t 310 \
    -c 48 \
    --resid \
    -m 0.001 \
    -p fluctmatch \
    -n 300 \
    --auto
```

---

## Run Modes

`run_fm` supports four execution behaviors:

- **fresh**  
  Default behavior when no mode flag is given. A new fluctuation-matching
  run is started from the initialized target state.

- **`--resume`**  
  Reloads the latest checkpointed run state, including parameter files,
  bond mask, and error history, and continues from the most conservative
  valid completed step.

- **`--restart`**  
  Rebuilds the fluctuation-matching state from the saved initialization
  files and restarts from cycle 1.

- **`--auto`**  
  Prioritizes resume when a valid checkpointed state exists, otherwise
  falls back to restart if initialization files are present, and otherwise
  starts fresh.

---

## What the Code Tracks During Fluctuation Matching

This implementation extends the original fluctuation-matching procedure
with bondwise convergence tracking, tail-convergence filtering, robust
fresh/restart/resume/auto execution modes, conservative checkpoint-based
recovery, and per-cycle timestamped progress logging.

### In Brief

1. Convergence is monitored bondwise rather than relying only on the
   overall root-mean-squared error.
2. At each cycle, every bond is evaluated individually against the
   fluctuation-difference tolerance.
3. The number and fraction of bonds that remain unconverged are tracked
   throughout the run using an active unconverged-bond mask.
4. When the remaining unconverged bonds fall below 0.3% of the total
   number of bonds in the system, a tail-convergence rule is applied.
5. In this tail region, a relative difference is evaluated as  
   `(fluct_diff - tol) / tol`.
6. Bonds with only a very small relative excess above tolerance are
   treated as effectively converged and removed from the active
   unconverged set.
7. Resume recovery is tolerant to minor checkpoint and `error.dat`
   step mismatches by continuing from the minimum consistent completed
   step rather than discarding recoverable progress.
8. Each completed fluctuation-matching cycle is logged with detailed
   timing information together with a wall-clock timestamp, making
   runtime progress easier to monitor on long HPC jobs.

---

## Typical Output Files Per Window

A completed window directory may contain files such as:

```text
data/1/
├── aa.xtc                  # split trajectory for this window
├── splittraj.log           # log file for trajectory splitting
├── cg.dcd                  # converted coarse-grained trajectory
├── cg.psf                  # converted coarse-grained structure file
├── fluctmatch.prm          # final fluctuation-matching parameters
├── fluctmatch.dist.prm     # final fluctuation-matching distance parameters
├── fluctmatch.xplor.psf    # final fluctuation-matching structure file
├── init.average.ic.        # target Internal coordinate averages
├── init.fluct.ic           # target Internal coordinate fluctuations
├── average.ic              # final average IC values across all cycles
├── fluct.ic                # final fluctuation IC values across all cycles
├── error.dat               # per-cycle error history
├── bond_mask.txt           # active unconverged bond mask across cycles
├── bond_convergence.txt    # bondwise convergence history across cycles
├── fm_checkpoint.json      # latest checkpointed state for recovery
└── charmm.log              # CHARMM log file for entire run

```

---

## Example Shell Loop Over All Windows

```bash
for win in ./data/*; do
    [ -d "$win" ] || continue

    fluctmatch convert \
        -s ./tests/data/system.tpr \
        -f "$win/aa.xtc" \
        -o "$win/" \
        -p cg \
        -m POLAR \ # if you have extra ligands add them as -M LIG (if CG for LIG is defined in the code base)
        -c 48 \
        --write

    fluctmatch convert \
        -s "$win/cg.psf" \
        -f "$win/cg.dcd" \
        -m ENM \
        --rmin 0.0 \
        --rmax 8.0 \
        -p "$win/fluctmatch" \
        -c 48

    fluctmatch run_fm \
        -s "$win/fluctmatch.xplor.psf" \
        -f "$win/cg.dcd" \
        -l "$win/charmm.log" \
        -o "$win/" \
        -e /path/to/charmm \
        -t 310 \
        -c 48 \
        --resid \
        -m 0.001 \
        -p fluctmatch \
        -n 300 \
        --auto
done
```

---

## Post-processing

After all window calculations are complete, parameter tables can be
consolidated across windows:

```bash
fluctmatch table \
    -d ./data \
    -o ./results \
    -p fluctmatch \
    -t Kb
```

If desired, the consolidated `Kb` table can then be converted back to
residue or bead labels for downstream analysis and network construction using `fluctmatch convert_table`.

