#!/usr/bin/env python3

"""Generate a GROMACS topology from a checked PDB file (gmx pdb2gmx).

Converted from nf-core/moleculardynamics modules/local/run_topology.nf.
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


def main():
    parser = argparse.ArgumentParser(description="Generate a GROMACS topology (gmx pdb2gmx)")
    parser.add_argument("--input", required=True, help="Checked PDB file")
    parser.add_argument("--output-gro", required=True, help="Output .gro structure file")
    parser.add_argument("--output-top", required=True, help="Output topology (.top) file")
    parser.add_argument("--output-itp", required=True, help="Output position restraint (.itp) file")
    parser.add_argument("--force-field", required=True, help="GROMACS force field (-ff)")
    parser.add_argument("--gmx-cmd", default="gmx", help="GROMACS command (default: gmx)")
    args = parser.parse_args()

    logger.info(f"Input: {args.input}")
    logger.info(f"Force field: {args.force_field}")

    for out in (args.output_gro, args.output_top, args.output_itp):
        out_dir = os.path.dirname(out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    cmd = [
        args.gmx_cmd, "pdb2gmx",
        "-f", args.input,
        "-o", args.output_gro,
        "-p", args.output_top,
        "-i", args.output_itp,
        "-ff", args.force_field,
        "-water", "spce",
        "-ignh",
    ]

    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    for out in (args.output_gro, args.output_top, args.output_itp):
        if not os.path.exists(out):
            logger.error(f"Expected output not found: {out}")
            sys.exit(1)

    logger.info(f"Outputs: {args.output_gro}, {args.output_top}, {args.output_itp}")


if __name__ == "__main__":
    main()
