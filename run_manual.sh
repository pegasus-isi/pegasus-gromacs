#!/bin/bash
#
# Manual test script for the GROMACS MD workflow.
# Runs each pipeline step locally (no Pegasus) against a small real test
# case — hen egg-white lysozyme (PDB 1AKI), the classic GROMACS
# protein-in-water tutorial system — to validate tool installation and
# argument wiring before submitting through Pegasus.
#
# Input data prep (download + .mdp files + samplesheet.csv) lives in
# prepare_test_data.sh, called below. Run that script on its own instead of
# this one when you just want input data ready for an actual Pegasus run.
#
# Usage:
#   ./run_manual.sh [--skip-download]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTPUT_DIR="${SCRIPT_DIR}/test_output"
GMX_CMD="${GMX_CMD:-gmx}"

SAMPLE="${SAMPLE:-1AKI}"
SAMPLE_DATA_DIR="${DATA_DIR}/${SAMPLE}"

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1"; }
log_step()    { echo ""; echo -e "${GREEN}========================================${NC}"; echo -e "${GREEN}STEP: $1${NC}"; echo -e "${GREEN}========================================${NC}"; }

mkdir -p "${OUTPUT_DIR}"

# ==============================================================
# Step 0: Prepare test data
# ==============================================================
"${SCRIPT_DIR}/prepare_test_data.sh" "$@"

PDB="${SAMPLE_DATA_DIR}/${SAMPLE}.pdb"

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
