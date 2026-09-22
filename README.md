# GROMACS MD Pegasus Workflow

A [Pegasus WMS](https://pegasus.isi.edu/) workflow that runs per-sample
GROMACS molecular dynamics simulations, converted from the
[nf-core/moleculardynamics](https://github.com/nf-core/moleculardynamics)
(nf-core/mdsimulations) Nextflow pipeline.

Each sample in the samplesheet runs an independent 10-step linear pipeline
(PDB cleanup through RMSD analysis and plotting). Samples have no cross-sample
dependencies, so Pegasus runs them all in parallel.

> **Note on the source pipeline:** nf-core/moleculardynamics is an early-stage
> scaffold with real bugs — several channel names in `workflows/mdsimulations.nf`
> are undefined, `RUN_SOLVATION`'s output tuple doesn't match
> `RUN_ENERGY_MINIMISATION`'s input arity, and `params.gmx_cmd` has no default.
> This conversion follows the *documented intent* of each module's script and
> description comment, not the broken Nextflow channel wiring. See
> "Conversion Notes" below for specifics.

## Pipeline Overview

```
structure.pdb ──> pdb_clean_and_check_missing_atoms ──> topology ──> solvation
                                                                         │
       ┌─────────────────────────────────────────────────────────────┘
       ▼
   energy_min ──> nvt_equilibration ──> npt_equilibration ──> production
                                                                  │
                                        ┌─────────────────────────┘
                                        ▼
                              post_processing ──> analysis_rmsd ──> analysis_plot ──> rmsd.png
                                                        │
                                                        └──> rmsd.xvg
```

| Step | Tool | Description |
|------|------|-------------|
| 1. pdb_clean_and_check_missing_atoms | `grep`-equivalent | Strip HETATM/CONECT records from the input PDB, then fail fast if it has missing atoms |
| 2. topology | `gmx pdb2gmx` | Generate topology (.gro/.top/.itp) from the PDB |
| 3. solvation | `gmx editconf/solvate/grompp/genion` | Build box, solvate, neutralize with ions |
| 4. energy_min | `gmx grompp/mdrun` | Energy minimization |
| 5. nvt_equilibration | `gmx grompp/mdrun` | NVT equilibration |
| 6. npt_equilibration | `gmx grompp/mdrun` | NPT equilibration |
| 7. production | `gmx grompp/mdrun/report-methods` | Production MD run + methods report |
| 8. post_processing | `gmx trjconv` | Remove periodicity (PBC) artifacts |
| 9. analysis_rmsd | `gmx rms` | RMSD of the trajectory (final output) |
| 10. analysis_plot | `matplotlib` | Plot the RMSD `.xvg` to a PNG (final output; not in the source pipeline) |

## Directory Structure

```
pegasus-gromacs/
├── workflow_generator.py         # Pegasus workflow generator
├── bin/
│   ├── pdb_clean_and_check_missing_atoms.py
│   ├── topology.py
│   ├── solvation.py
│   ├── energy_min.py
│   ├── nvt_equilibration.py
│   ├── npt_equilibration.py
│   ├── production.py
│   ├── post_processing.py
│   ├── analysis_rmsd.py
│   └── analysis_plot.py
├── Apptainer/
│   └── Gromacs_Container.def     # Container definition (micromamba + GROMACS + matplotlib)
├── data/
│   ├── samplesheet.csv           # Samplesheet (written by run_manual.sh)
│   └── test/                     # Test input data (created by run_manual.sh)
│       └── <sample>/             # One subfolder per sample, e.g. 1AKI/
├── run_manual.sh                 # Local smoke test, no Pegasus required
└── README.md
```

## Prerequisites

- [Pegasus WMS](https://pegasus.isi.edu/) >= 5.0
- [HTCondor](https://htcondor.org/) >= 10.2
- Python 3.8+ with `pyyaml` installed (needed by the Pegasus Python API)
- [Apptainer](https://apptainer.org/) (to build the GROMACS container)

## Setup

### 1. Build the Apptainer Container

```bash
cd pegasus-gromacs
apptainer build Apptainer/Gromacs_Container.sif Apptainer/Gromacs_Container.def
```

The `image=` path in `workflow_generator.py`'s `Container()` definition
already points at `Apptainer/Gromacs_Container.sif` relative to the workflow
directory — no edit needed if you build it there.

### 2. Prepare Input Data

Each sample needs: a structure file (PDB), and four `.mdp` files (energy
minimization, NVT equilibration, NPT equilibration, production MD). Keep each
sample's files in their own subfolder, named after the sample — not a shared
folder — so multiple samples' PDB/mdp files never collide. List them in a
samplesheet CSV:

```csv
sample,structure,em_mdp,nvt_mdp,npt_mdp,md_mdp,force_field,box_type,distance_to_box
1AKI,data/test/1AKI/1AKI.pdb,data/test/1AKI/em.mdp,data/test/1AKI/nvt.mdp,data/test/1AKI/npt.mdp,data/test/1AKI/md.mdp,charmm27,cubic,1.0
```

- `force_field`: one of `charmm27`, `charmm36`, `amber`, `amber99sb`, `amber14sb`
- `box_type` (optional, default `cubic`): `cubic`, `triclinic`, or `dodecahedron`
- `distance_to_box` (optional, default `1.0`): distance to box edge in nm, `0.0`–`5.0`

To generate a small real test case (hen egg-white lysozyme, PDB 1AKI) with
placeholder `.mdp` files under `data/test/1AKI/`, and write a matching
`data/samplesheet.csv`, run:

```bash
./run_manual.sh
```

## Usage

### Generate Workflow

```bash
python3 workflow_generator.py --samplesheet data/samplesheet.csv --output workflow.yml
```

### CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `--samplesheet` | (required) | CSV samplesheet (see format above) |
| `--gmx-cmd` | `gmx` | GROMACS command to invoke (e.g. `gmx_mpi`) |
| `-e`, `--execution-site-name` | `condorpool` | HTCondor execution site name |
| `-s`, `--skip-sites-catalog` | false | Skip site catalog creation |
| `-o`, `--output` | `workflow.yml` | Output workflow file |

### Plan and Submit

```bash
pegasus-plan --dir submit -s condorpool -o local workflow.yml
```

### Monitor

```bash
pegasus-status <run-directory>
pegasus-statistics <run-directory>
```

## Outputs

Final outputs staged to the `output/` directory, per sample:

| Output | Description |
|--------|-------------|
| `{sample}_md.gro`, `{sample}_md.tpr`, `{sample}_md.xtc` | Production MD trajectory + structure |
| `{sample}_MD_REPORT.out` | Simulation methods report (`gmx report-methods`; `-o` requires a `.out` suffix) |
| `{sample}_rmsd.xvg` | RMSD of the trajectory — final analysis output |
| `{sample}_rmsd.png` | RMSD vs. time plot (not in the source pipeline) |

## Resource Requirements

| Step | Memory | Cores |
|------|--------|-------|
| pdb_clean_and_check_missing_atoms | 1 GB | 1 |
| topology | 2 GB | 1 |
| solvation | 2 GB | 1 |
| energy_min | 4 GB | 4 |
| nvt_equilibration | 4 GB | 4 |
| npt_equilibration | 4 GB | 4 |
| production | 8 GB | 4 |
| post_processing | 2 GB | 1 |
| analysis_rmsd | 2 GB | 1 |
| analysis_plot | 1 GB | 1 |

## Conversion Notes (Nextflow → Pegasus)

| nf-core/moleculardynamics | Pegasus |
|----------------------------|---------|
| `PRE_POS_CLEAN_PDB` + `PRE_POS_CHECK_MISSING_ATOMS` processes | Combined into a single `pdb_clean_and_check_missing_atoms` Transformation + Job |
| `RUN_TOPOLOGY` process | `topology` Transformation + Job |
| `RUN_SOLVATION` process | `solvation` Transformation + Job |
| `RUN_ENERGY_MINIMISATION` process | `energy_min` Transformation + Job |
| `RUN_NVT_EQUILIBRATION` process | `nvt_equilibration` Transformation + Job |
| `RUN_NPT_EQUILIBRATION` process | `npt_equilibration` Transformation + Job |
| `RUN_PRODUCTION` process | `production` Transformation + Job |
| `POST_PROCESSING` process | `post_processing` Transformation + Job |
| `ANALYSIS_RMSD` process | `analysis_rmsd` Transformation + Job |
| *(none — added on top of the conversion)* | `analysis_plot` Transformation + Job: plots `rmsd.xvg` to `rmsd.png` with matplotlib |
| `--input` samplesheet (nf-schema) | `--samplesheet` CSV, parsed with `csv.DictReader` |
| Per-sample channel tuples | Python `File()` objects shared between jobs |
| `params.gmx_cmd` (no default upstream) | `--gmx-cmd` CLI arg, default `gmx` |
| Implicit `topol.top` mutated in place by `genion` | Explicit `{sample}_topol_solv.top`, a distinct `File` produced by `solvation` and passed read-only to later steps |
| `${mdp.simpleName}` as GROMACS `-deffnm` | Explicit `--output-gro` per step; wrapper derives `-deffnm` from it, guaranteeing unique per-sample output names |

Notable deviations from the literal (broken) source wiring:
- `RUN_SOLVATION`'s and `RUN_ENERGY_MINIMISATION`'s channel tuples don't line
  up in the source (different arity/order) and would not actually run in
  Nextflow as checked in. The Pegasus `solvation` → `energy_min` job pair
  wires the same *data* (solvated+ionized `.gro`, updated topology, itp) via
  explicit named arguments instead.
- `workflows/mdsimulations.nf` references undefined channels
  (`ch_postprocessed_gro`, `ch_postprocessed_tpr`, `ch_cleaned_pdb`). The
  Pegasus DAG instead follows each module's `script:` block and header
  comment literally: `post_processing` consumes the production `.tpr`/`.xtc`
  (matching `POST_PROCESSING`'s actual `trjconv` command), and
  `analysis_rmsd` consumes the production `.gro` and the PBC-corrected `.xtc`
  (matching `ANALYSIS_RMSD`'s actual `gmx rms -s <gro> -f <xtc>` command).
