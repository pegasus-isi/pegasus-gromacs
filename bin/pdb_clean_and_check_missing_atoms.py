#!/usr/bin/env python3

"""Strip HETATM/CONECT records from a PDB file, then fail if it has missing atoms.

Combines the two source pipeline steps into a single job:
clean_pdb (nf-core/moleculardynamics modules/local/pre_pos_clean_pdb.nf) and
check_missing_atoms (modules/local/pre_pos_check_missing_atoms.nf).
"""

import argparse
import logging
import os
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Strip HETATM/CONECT records from a PDB file and check for missing atoms"
    )
    parser.add_argument("--input", required=True, help="Input PDB file")
    parser.add_argument("--output", required=True, help="Cleaned, checked PDB file")
    args = parser.parse_args()

    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.input) as fh:
        lines = fh.readlines()

    cleaned = [
        line for line in lines
        if not line.startswith("HETATM") and not line.startswith("CONECT")
    ]
    logger.info(f"Removed {len(lines) - len(cleaned)} HETATM/CONECT line(s)")

    missing_lines = [line.rstrip() for line in cleaned if "MISSING" in line]
    if missing_lines:
        logger.error(f"Missing atoms found in {args.input}")
        for line in missing_lines:
            print(line, file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w") as fh:
        fh.writelines(cleaned)

    if not os.path.exists(args.output):
        logger.error(f"Expected output not found: {args.output}")
        sys.exit(1)

    logger.info(f"No missing atoms found. Output: {args.output}")


if __name__ == "__main__":
    main()
