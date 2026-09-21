#!/usr/bin/env python3

"""Remove periodicity artifacts from the production trajectory (gmx trjconv).

Converted from nf-core/moleculardynamics modules/local/post_processing.nf.
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
    parser = argparse.ArgumentParser(description="Remove PBC artifacts from an MD trajectory")
    parser.add_argument("--input-tpr", required=True, help="Production .tpr file")
    parser.add_argument("--input-xtc", required=True, help="Production trajectory (.xtc)")
    parser.add_argument("--output", required=True, help="Output PBC-corrected trajectory (.xtc)")
    parser.add_argument("--gmx-cmd", default="gmx", help="GROMACS command (default: gmx)")
    args = parser.parse_args()

    logger.info(f"Input tpr: {args.input_tpr}, xtc: {args.input_xtc}")

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    cmd = [
        args.gmx_cmd, "trjconv",
        "-s", args.input_tpr,
        "-f", args.input_xtc,
        "-o", args.output,
        "-pbc", "mol",
        "-center",
    ]

    logger.info(f"Running: {' '.join(cmd)}")
    # Select group 1 (protein) to center, group 0 (system) to output — matches
    # the source pipeline's `printf "1\n0\n"`.
    result = subprocess.run(cmd, input="1\n0\n", capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    if not os.path.exists(args.output):
        logger.error(f"Expected output not found: {args.output}")
        sys.exit(1)

    logger.info(f"Post-processing completed! Output: {args.output}")


if __name__ == "__main__":
    main()
