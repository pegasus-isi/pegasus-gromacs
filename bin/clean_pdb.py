#!/usr/bin/env python3

"""Strip HETATM and CONECT records from a PDB file.

Converted from nf-core/moleculardynamics modules/local/pre_pos_clean_pdb.nf.
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
    parser = argparse.ArgumentParser(description="Remove HETATM/CONECT records from a PDB file")
    parser.add_argument("--input", required=True, help="Input PDB file")
    parser.add_argument("--output", required=True, help="Cleaned PDB file")
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

    with open(args.output, "w") as fh:
        fh.writelines(cleaned)

    if not os.path.exists(args.output):
        logger.error(f"Expected output not found: {args.output}")
        sys.exit(1)

    logger.info(f"Removed {len(lines) - len(cleaned)} HETATM/CONECT line(s)")
    logger.info(f"Output: {args.output}")


if __name__ == "__main__":
    main()
