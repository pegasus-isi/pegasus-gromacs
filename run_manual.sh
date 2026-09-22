#!/bin/bash
#
# Manual test script for the GROMACS MD workflow.
# Runs each pipeline step locally (no Pegasus) against a small real test
# case — hen egg-white lysozyme (PDB 1AKI), the classic GROMACS
# protein-in-water tutorial system — to validate tool installation and
# argument wiring before submitting through Pegasus. Also (re)writes
# data/samplesheet.csv so the same downloaded/generated data can be fed
# straight into workflow_generator.py for an actual Pegasus run.
#
# Each sample's input data lives in its own subfolder, named after the
# sample (data/test/<sample>/), not a shared "test" folder — so multiple
# samples' PDB/mdp files never collide.
#
# The .mdp files generated below use deliberately tiny step counts so the
# whole chain finishes in well under a minute. They are NOT suitable for a
# real production MD run — replace them with properly validated .mdp files
# (e.g. from the GROMACS tutorials) for real simulations.
#
# Usage:
#   ./run_manual.sh [--skip-download]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DATA_DIR="${SCRIPT_DIR}/data/"
OUTPUT_DIR="${SCRIPT_DIR}/test_output"
SAMPLESHEET="${SCRIPT_DIR}/data/samplesheet.csv"
GMX_CMD="${GMX_CMD:-gmx}"

SAMPLE="1AKI"
SAMPLE_DATA_DIR="${TEST_DATA_DIR}/${SAMPLE}"

SKIP_DOWNLOAD=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-download) SKIP_DOWNLOAD=true; shift ;;
        *) echo "Unknown argument: $1"; echo "Usage: $0 [--skip-download]"; exit 1 ;;
    esac
done

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1"; }
log_step()    { echo ""; echo -e "${GREEN}========================================${NC}"; echo -e "${GREEN}STEP: $1${NC}"; echo -e "${GREEN}========================================${NC}"; }

mkdir -p "${SAMPLE_DATA_DIR}" "${OUTPUT_DIR}"

# ==============================================================
# Step 0: Prepare test data
# ==============================================================
log_step "Preparing test data"

PDB="${SAMPLE_DATA_DIR}/${SAMPLE}.pdb"

if [ "$SKIP_DOWNLOAD" = false ]; then
    if [ ! -f "${PDB}" ]; then
        log_info "Downloading ${SAMPLE}.pdb (hen egg-white lysozyme) from RCSB..."
        curl -sL -o "${PDB}" "https://files.rcsb.org/download/${SAMPLE}.pdb"
    fi
else
    log_info "Skipping download (--skip-download)"
fi

# Minimal .mdp files — short nsteps for a fast smoke test only.
cat > "${SAMPLE_DATA_DIR}/em.mdp" <<'EOF'
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF

cat > "${SAMPLE_DATA_DIR}/nvt.mdp" <<'EOF'
integrator  = md
nsteps      = 50
dt          = 0.002
continuation = no
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
pbc         = xyz
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1
EOF

cat > "${SAMPLE_DATA_DIR}/npt.mdp" <<'EOF'
integrator  = md
nsteps      = 50
dt          = 0.002
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
pcoupl      = C-rescale
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
gen_vel     = no
EOF

cat > "${SAMPLE_DATA_DIR}/md.mdp" <<'EOF'
integrator  = md
nsteps      = 100
dt          = 0.002
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
pcoupl      = C-rescale
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
nstxout-compressed = 10
EOF

log_success "Test data ready in ${SAMPLE_DATA_DIR}"

mkdir -p "$(dirname "${SAMPLESHEET}")"
cat > "${SAMPLESHEET}" <<EOF
sample,structure,em_mdp,nvt_mdp,npt_mdp,md_mdp,force_field,box_type,distance_to_box
${SAMPLE},data/test/${SAMPLE}/${SAMPLE}.pdb,data/test/${SAMPLE}/em.mdp,data/test/${SAMPLE}/nvt.mdp,data/test/${SAMPLE}/npt.mdp,data/test/${SAMPLE}/md.mdp,charmm27,cubic,1.0
EOF
log_success "Samplesheet written to ${SAMPLESHEET}"

# ==============================================================
# Step 1: pdb_clean_and_check_missing_atoms
# ==============================================================
log_step "1. pdb_clean_and_check_missing_atoms"
python3 "${SCRIPT_DIR}/bin/pdb_clean_and_check_missing_atoms.py" \
    --input "${PDB}" \
    --output "${OUTPUT_DIR}/${SAMPLE}_checked.pdb"

# ==============================================================
# Step 2: topology (gmx pdb2gmx)
# ==============================================================
log_step "2. topology"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/topology.py" \
        --input "${SAMPLE}_checked.pdb" \
        --output-gro "${SAMPLE}_topology.gro" \
        --output-top "${SAMPLE}_topol.top" \
        --output-itp "${SAMPLE}_posre.itp" \
        --force-field charmm27 \
        --gmx-cmd "${GMX_CMD}"
)

# ==============================================================
# Step 3: solvation
# ==============================================================
log_step "3. solvation"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/solvation.py" \
        --input-gro "${SAMPLE}_topology.gro" \
        --input-top "${SAMPLE}_topol.top" \
        --output-gro "${SAMPLE}_box_solv_ions.gro" \
        --output-top "${SAMPLE}_topol_solv.top" \
        --box-type cubic \
        --distance-to-box 1.0 \
        --gmx-cmd "${GMX_CMD}"
)

# ==============================================================
# Step 4: energy_min
# ==============================================================
log_step "4. energy_min"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/energy_min.py" \
        --input-gro "${SAMPLE}_box_solv_ions.gro" \
        --input-top "${SAMPLE}_topol_solv.top" \
        --mdp "${SAMPLE_DATA_DIR}/em.mdp" \
        --output-gro "${SAMPLE}_em.gro" \
        --gmx-cmd "${GMX_CMD}"
)

# ==============================================================
# Step 5: nvt_equilibration
# ==============================================================
log_step "5. nvt_equilibration"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/nvt_equilibration.py" \
        --input-gro "${SAMPLE}_em.gro" \
        --input-top "${SAMPLE}_topol_solv.top" \
        --mdp "${SAMPLE_DATA_DIR}/nvt.mdp" \
        --output-gro "${SAMPLE}_nvt.gro" \
        --gmx-cmd "${GMX_CMD}"
)

# ==============================================================
# Step 6: npt_equilibration
# ==============================================================
log_step "6. npt_equilibration"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/npt_equilibration.py" \
        --input-gro "${SAMPLE}_nvt.gro" \
        --input-top "${SAMPLE}_topol_solv.top" \
        --mdp "${SAMPLE_DATA_DIR}/npt.mdp" \
        --output-gro "${SAMPLE}_npt.gro" \
        --gmx-cmd "${GMX_CMD}"
)

# ==============================================================
# Step 7: production
# ==============================================================
log_step "7. production"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/production.py" \
        --input-gro "${SAMPLE}_npt.gro" \
        --input-top "${SAMPLE}_topol_solv.top" \
        --mdp "${SAMPLE_DATA_DIR}/md.mdp" \
        --output-gro "${SAMPLE}_md.gro" \
        --output-tpr "${SAMPLE}_md.tpr" \
        --output-xtc "${SAMPLE}_md.xtc" \
        --output-report "${SAMPLE}_MD_REPORT.out" \
        --gmx-cmd "${GMX_CMD}"
)

# ==============================================================
# Step 8: post_processing
# ==============================================================
log_step "8. post_processing"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/post_processing.py" \
        --input-tpr "${SAMPLE}_md.tpr" \
        --input-xtc "${SAMPLE}_md.xtc" \
        --output "${SAMPLE}_noPBC.xtc" \
        --gmx-cmd "${GMX_CMD}"
)

# ==============================================================
# Step 9: analysis_rmsd
# ==============================================================
log_step "9. analysis_rmsd"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/analysis_rmsd.py" \
        --input-gro "${SAMPLE}_md.gro" \
        --input-xtc "${SAMPLE}_noPBC.xtc" \
        --output "${SAMPLE}_rmsd.xvg" \
        --gmx-cmd "${GMX_CMD}"
)

# ==============================================================
# Step 10: analysis_plot (RMSD .xvg -> PNG; not in the source pipeline)
# ==============================================================
log_step "10. analysis_plot"
(
    cd "${OUTPUT_DIR}" && python3 "${SCRIPT_DIR}/bin/analysis_plot.py" \
        --input "${SAMPLE}_rmsd.xvg" \
        --output "${SAMPLE}_rmsd.png" \
        --title "${SAMPLE} RMSD vs. Time"
)

echo ""
echo "=============================================="
echo "  TEST COMPLETED SUCCESSFULLY!"
echo "=============================================="
echo ""
ls -lh "${OUTPUT_DIR}"
echo ""
log_success "All steps passed! Ready to run with Pegasus."
