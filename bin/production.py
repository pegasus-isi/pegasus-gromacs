#!/usr/bin/env python3

"""Run the production MD simulation and generate a methods report.

Runs gmx grompp, gmx mdrun, and gmx report-methods.
Converted from nf-core/moleculardynamics modules/local/run_production.nf.
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
    parser = argparse.ArgumentParser(description="Run GROMACS production MD")
    parser.add_argument("--input-gro", required=True, help="NPT-equilibrated .gro structure")
    parser.add_argument("--input-top", required=True, help="Topology (.top) file")
    parser.add_argument("--mdp", required=True, help="Production MD .mdp file")
    parser.add_argument("--output-gro", required=True, help="Output production .gro file")
    parser.add_argument("--output-tpr", required=True, help="Output production .tpr file")
    parser.add_argument("--output-xtc", required=True, help="Output production trajectory (.xtc)")
    parser.add_argument("--output-report", required=True, help="Output MD_REPORT file")
    parser.add_argument("--gmx-cmd", default="gmx", help="GROMACS command (default: gmx)")
    args = parser.parse_args()

    logger.info(f"Input: {args.input_gro}")

    for out in (args.output_gro, args.output_tpr, args.output_xtc, args.output_report):
        out_dir = os.path.dirname(out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    deffnm = os.path.splitext(os.path.basename(args.output_gro))[0]

    run([args.gmx_cmd, "grompp", "-f", args.mdp, "-c", args.input_gro,
         "-p", args.input_top, "-o", f"{deffnm}.tpr"])
    run([args.gmx_cmd, "mdrun", "-v", "-deffnm", deffnm])
    run([args.gmx_cmd, "report-methods", "-s", f"{deffnm}.tpr", "-o", args.output_report])

    renames = {
        f"{deffnm}.gro": args.output_gro,
        f"{deffnm}.tpr": args.output_tpr,
        f"{deffnm}.xtc": args.output_xtc,
    }
    for produced, declared in renames.items():
        if produced != declared:
            os.replace(produced, declared)

    for out in (args.output_gro, args.output_tpr, args.output_xtc, args.output_report):
        if not os.path.exists(out):
            logger.error(f"Expected output not found: {out}")
            sys.exit(1)

    logger.info(f"Simulation completed! Outputs: {args.output_gro}, {args.output_tpr}, "
                f"{args.output_xtc}, {args.output_report}")


if __name__ == "__main__":
    main()
