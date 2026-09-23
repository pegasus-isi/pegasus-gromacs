#!/bin/bash
#
# Prepares real, small GROMACS MD input data — hen egg-white lysozyme
# (PDB 1AKI), the classic GROMACS protein-in-water tutorial system — and
# writes data/samplesheet.csv pointing at it.
#
# Run this on its own to get input data ready for an actual Pegasus run
# (workflow_generator.py --samplesheet data/samplesheet.csv ...). It's also
# called by run_manual.sh before its local, no-Pegasus smoke test.
#
# Each sample's input data lives in its own subfolder, named after the
# sample (data/<sample>/), not a shared folder — so multiple samples'
# PDB/mdp files never collide.
#
# The .mdp files written below use deliberately tiny step counts so a local
# smoke test finishes in well under a minute. They are NOT suitable for a
# real production MD run — replace them with properly validated .mdp files
# (e.g. from the GROMACS tutorials) for real simulations.
#
# Usage:
#   ./prepare_test_data.sh [--skip-download]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
SAMPLESHEET="${DATA_DIR}/samplesheet.csv"

SAMPLE="${SAMPLE:-1AKI}"
SAMPLE_DATA_DIR="${DATA_DIR}/${SAMPLE}"

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
log_step()    { echo ""; echo -e "${GREEN}========================================${NC}"; echo -e "${GREEN}STEP: $1${NC}"; echo -e "${GREEN}========================================${NC}"; }

mkdir -p "${SAMPLE_DATA_DIR}"

log_step "Preparing test data for sample: ${SAMPLE}"

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
${SAMPLE},data/${SAMPLE}/${SAMPLE}.pdb,data/${SAMPLE}/em.mdp,data/${SAMPLE}/nvt.mdp,data/${SAMPLE}/npt.mdp,data/${SAMPLE}/md.mdp,charmm27,cubic,1.0
EOF
log_success "Samplesheet written to ${SAMPLESHEET}"

echo ""
log_info "To run the actual workflow with Pegasus:"
log_info "  python3 workflow_generator.py --samplesheet ${SAMPLESHEET#${SCRIPT_DIR}/} --output workflow.yml"
