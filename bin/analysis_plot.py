#!/usr/bin/env python3

"""Plot an RMSD .xvg time series (gmx rms output) to a PNG image.

Not part of the source nf-core/moleculardynamics pipeline — added on top of
the conversion so the RMSD result can be inspected without a separate
plotting tool (xmgrace, etc.).
"""

import argparse
import logging
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def read_xvg(path):
    """Parse a GROMACS .xvg file, ignoring comment/metadata lines."""
    xs, ys = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
    return xs, ys


def main():
    parser = argparse.ArgumentParser(description="Plot an RMSD .xvg file to PNG")
    parser.add_argument("--input", required=True, help="RMSD data file (.xvg)")
    parser.add_argument("--output", required=True, help="Output plot image (.png)")
    parser.add_argument("--title", default="RMSD vs. Time", help="Plot title")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)

    xs, ys = read_xvg(args.input)
    if not xs:
        logger.error(f"No data points parsed from {args.input}")
        sys.exit(1)

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xs, ys, linewidth=1.2)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (nm)")
    ax.set_title(args.title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)

    if not os.path.exists(args.output):
        logger.error(f"Expected output not found: {args.output}")
        sys.exit(1)

    logger.info(f"RMSD plot written to {args.output}")


if __name__ == "__main__":
    main()
