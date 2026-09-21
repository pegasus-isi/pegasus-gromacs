#!/usr/bin/env python3

"""Calculate the RMSD of the protein along the MD trajectory (gmx rms).

Converted from nf-core/moleculardynamics modules/local/analysis_rmsd.nf.
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
    parser = argparse.ArgumentParser(description="Calculate RMSD of an MD trajectory")
    parser.add_argument("--input-gro", required=True, help="Structure reference (.gro)")
    parser.add_argument("--input-xtc", required=True, help="PBC-corrected trajectory (.xtc)")
    parser.add_argument("--output", required=True, help="Output RMSD data file (.xvg)")
    parser.add_argument("--gmx-cmd", default="gmx", help="GROMACS command (default: gmx)")
    args = parser.parse_args()

    logger.info(f"Input gro: {args.input_gro}, xtc: {args.input_xtc}")

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    cmd = [
        args.gmx_cmd, "rms",
        "-s", args.input_gro,
        "-f", args.input_xtc,
        "-o", args.output,
        "-tu", "ns",
    ]

    logger.info(f"Running: {' '.join(cmd)}")
    # Select group 3 (C-alpha) as the fit group, group 1 (protein) as the
    # analysis group — matches the source pipeline's `printf "3\n1\n"`.
    result = subprocess.run(cmd, input="3\n1\n", capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    if not os.path.exists(args.output):
        logger.error(f"Expected output not found: {args.output}")
        sys.exit(1)

    logger.info(f"RMSD analysis completed! Output: {args.output}")


if __name__ == "__main__":
    main()
