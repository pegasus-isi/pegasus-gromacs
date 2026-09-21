#!/usr/bin/env python3

"""Run energy minimization (gmx grompp + gmx mdrun).

Converted from nf-core/moleculardynamics modules/local/run_energy_min.nf.
"""

import argparse
import logging
import os
import subprocess
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def run(cmd):
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)
    return result


def main():
    parser = argparse.ArgumentParser(description="Run GROMACS energy minimization")
    parser.add_argument("--input-gro", required=True, help="Solvated+ionized .gro structure")
    parser.add_argument("--input-top", required=True, help="Topology (.top) file")
    parser.add_argument("--mdp", required=True, help="Energy minimization .mdp file")
    parser.add_argument("--output-gro", required=True, help="Output minimized .gro file")
    parser.add_argument("--gmx-cmd", default="gmx", help="GROMACS command (default: gmx)")
    args = parser.parse_args()

    logger.info(f"Input: {args.input_gro}")

    out_dir = os.path.dirname(args.output_gro)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    deffnm = os.path.splitext(os.path.basename(args.output_gro))[0]
    tpr = f"{deffnm}.tpr"

    run([args.gmx_cmd, "grompp", "-f", args.mdp, "-c", args.input_gro,
         "-p", args.input_top, "-o", tpr])
    run([args.gmx_cmd, "mdrun", "-v", "-deffnm", deffnm])

    produced_gro = f"{deffnm}.gro"
    if produced_gro != args.output_gro:
        os.replace(produced_gro, args.output_gro)

    if not os.path.exists(args.output_gro):
        logger.error(f"Expected output not found: {args.output_gro}")
        sys.exit(1)

    logger.info("Energy minimization completed")
    logger.info(f"Output: {args.output_gro}")


if __name__ == "__main__":
    main()
