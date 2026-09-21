#!/usr/bin/env python3

"""Fail the workflow if a cleaned PDB file has missing atoms.

Converted from nf-core/moleculardynamics
modules/local/pre_pos_check_missing_atoms.nf.
"""

import argparse
import logging
import os
import shutil
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Check for missing atoms in a PDB file")
    parser.add_argument("--input", required=True, help="Cleaned PDB file")
    parser.add_argument("--output", required=True, help="Checked PDB file (copy of input)")
    args = parser.parse_args()

    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    missing_lines = []
    with open(args.input) as fh:
        for line in fh:
            if "MISSING" in line:
                missing_lines.append(line.rstrip())

    if missing_lines:
        logger.error(f"Missing atoms found in {args.input}")
        for line in missing_lines:
            print(line, file=sys.stderr)
        sys.exit(1)

    shutil.copy(args.input, args.output)
    logger.info(f"No missing atoms found in {args.input}")

    if not os.path.exists(args.output):
        logger.error(f"Expected output not found: {args.output}")
        sys.exit(1)

    logger.info(f"Output: {args.output}")


if __name__ == "__main__":
    main()
