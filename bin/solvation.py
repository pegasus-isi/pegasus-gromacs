#!/usr/bin/env python3

"""Define the simulation box, solvate, and neutralize with ions.

Runs gmx editconf, gmx solvate, gmx grompp, and gmx genion.
Converted from nf-core/moleculardynamics modules/local/run_solvation.nf.

genion mutates the topology file in place (updates molecule counts), so this
wrapper copies the input topology to a local working copy, runs the GROMACS
chain against it, then copies the mutated result to the declared output path.
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def run(cmd, gmx_cmd, **kwargs):
    logger.info(f"Running: {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)
    return result


def main():
    parser = argparse.ArgumentParser(description="Build box, solvate, and add ions")
    parser.add_argument("--input-gro", required=True, help="Input .gro structure (from topology step)")
    parser.add_argument("--input-top", required=True, help="Input topology (.top) file")
    parser.add_argument("--output-gro", required=True, help="Output solvated+ionized .gro file")
    parser.add_argument("--output-top", required=True, help="Output updated topology (.top) file")
    parser.add_argument("--box-type", default="cubic", help="Box type (default: cubic)")
    parser.add_argument("--distance-to-box", default="1.0", help="Distance to box edge in nm")
    parser.add_argument("--gmx-cmd", default="gmx", help="GROMACS command (default: gmx)")
    args = parser.parse_args()

    logger.info(f"Input gro: {args.input_gro}")
    logger.info(f"Box type: {args.box_type}, distance: {args.distance_to_box} nm")

    for out in (args.output_gro, args.output_top):
        out_dir = os.path.dirname(out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    # Local working copy of the topology — genion/grompp mutate it in place.
    working_top = "working_topol.top"
    shutil.copy(args.input_top, working_top)

    sample_prefix = os.path.splitext(os.path.basename(args.output_gro))[0]
    box_gro = f"{sample_prefix}_box.gro"
    box_solv_gro = f"{sample_prefix}_box_solv.gro"
    ions_mdp = "ions.mdp"
    ions_tpr = f"{sample_prefix}_ions.tpr"

    run(
        [args.gmx_cmd, "editconf", "-f", args.input_gro, "-o", box_gro,
         "-c", "-d", args.distance_to_box, "-bt", args.box_type],
        args.gmx_cmd,
    )
    run(
        [args.gmx_cmd, "solvate", "-cp", box_gro, "-cs", "spc216.gro",
         "-o", box_solv_gro, "-p", working_top],
        args.gmx_cmd,
    )

    with open(ions_mdp, "w") as fh:
        pass  # empty placeholder .mdp, matching the source pipeline's `touch ions.mdp`

    run(
        [args.gmx_cmd, "grompp", "-f", ions_mdp, "-c", box_solv_gro,
         "-p", working_top, "-o", ions_tpr],
        args.gmx_cmd,
    )
    run(
        [args.gmx_cmd, "genion", "-s", ions_tpr, "-o", args.output_gro,
         "-p", working_top, "-neutral", "-conc", "0.15", "-pname", "NA", "-nname", "CL"],
        args.gmx_cmd,
        input="SOL\n",
    )

    shutil.copy(working_top, args.output_top)

    for out in (args.output_gro, args.output_top):
        if not os.path.exists(out):
            logger.error(f"Expected output not found: {out}")
            sys.exit(1)

    logger.info(f"Outputs: {args.output_gro}, {args.output_top}")


if __name__ == "__main__":
    main()
