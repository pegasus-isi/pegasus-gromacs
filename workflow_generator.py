#!/usr/bin/env python3

"""
Pegasus workflow generator for GROMACS molecular dynamics simulations.

Converted from the nf-core/moleculardynamics (nf-core/mdsimulations) Nextflow
pipeline: https://github.com/nf-core/moleculardynamics

Each sample in the samplesheet runs an independent, linear 10-step GROMACS MD
pipeline (structure cleanup through RMSD analysis). Samples are fully
independent, so they all run in parallel.

Pipeline steps (per sample):
1. pdb_clean_and_check_missing_atoms - strip HETATM/CONECT records from the
                           input PDB, then fail fast if it has missing atoms
2. topology              - gmx pdb2gmx: PDB -> topology (.gro/.top/.itp)
3. solvation             - gmx editconf/solvate/grompp/genion: box + solvent + ions
4. energy_min            - gmx grompp/mdrun: energy minimization
5. nvt_equilibration     - gmx grompp/mdrun: NVT equilibration
6. npt_equilibration     - gmx grompp/mdrun: NPT equilibration
7. production            - gmx grompp/mdrun/report-methods: production MD run
8. post_processing       - gmx trjconv: remove periodicity artifacts
9. analysis_rmsd         - gmx rms: RMSD of the trajectory
10. analysis_plot        - matplotlib: plot the RMSD .xvg to a PNG
                           (not in the source pipeline; added for visualization)

Usage:
    ./workflow_generator.py --samplesheet samplesheet.csv --output workflow.yml
"""

import argparse
import csv
import logging
import os
import sys
from pathlib import Path

from Pegasus.api import *

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


VALID_FORCE_FIELDS = {"charmm27", "charmm36", "amber", "amber99sb", "amber14sb"}
VALID_BOX_TYPES = {"cubic", "triclinic", "dodecahedron"}

REQUIRED_COLUMNS = [
    "sample",
    "structure",
    "em_mdp",
    "nvt_mdp",
    "npt_mdp",
    "md_mdp",
    "force_field",
]

# Per-tool resource configuration for the transformation catalog.
# mdrun steps (energy_min, nvt/npt equilibration, production) are the
# expensive ones; the rest are cheap bookkeeping/analysis steps.
TOOL_CONFIGS = {
    "pdb_clean_and_check_missing_atoms": {"memory": "1 GB", "cores": 1},
    "topology": {"memory": "2 GB", "cores": 1},
    "solvation": {"memory": "2 GB", "cores": 1},
    "energy_min": {"memory": "4 GB", "cores": 4},
    "nvt_equilibration": {"memory": "4 GB", "cores": 4},
    "npt_equilibration": {"memory": "4 GB", "cores": 4},
    "production": {"memory": "8 GB", "cores": 4},
    "post_processing": {"memory": "2 GB", "cores": 1},
    "analysis_rmsd": {"memory": "2 GB", "cores": 1},
    "analysis_plot": {"memory": "1 GB", "cores": 1},
}


def parse_samplesheet(path):
    """Parse the samplesheet CSV into a list of validated sample dicts."""
    samples = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing_cols = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing_cols:
            print(
                f"Error: samplesheet is missing required column(s): {missing_cols}",
                file=sys.stderr,
            )
            sys.exit(1)

        for row in reader:
            sample = {
                "sample": row["sample"].strip(),
                "structure": row["structure"].strip(),
                "em_mdp": row["em_mdp"].strip(),
                "nvt_mdp": row["nvt_mdp"].strip(),
                "npt_mdp": row["npt_mdp"].strip(),
                "md_mdp": row["md_mdp"].strip(),
                "force_field": row["force_field"].strip(),
                "box_type": (row.get("box_type") or "cubic").strip() or "cubic",
                "distance_to_box": (row.get("distance_to_box") or "1.0").strip()
                or "1.0",
            }

            if not sample["sample"]:
                print("Error: samplesheet row has an empty 'sample' value", file=sys.stderr)
                sys.exit(1)
            if sample["force_field"] not in VALID_FORCE_FIELDS:
                print(
                    f"Error: sample {sample['sample']!r} has invalid force_field "
                    f"{sample['force_field']!r} (must be one of {sorted(VALID_FORCE_FIELDS)})",
                    file=sys.stderr,
                )
                sys.exit(1)
            if sample["box_type"] not in VALID_BOX_TYPES:
                print(
                    f"Error: sample {sample['sample']!r} has invalid box_type "
                    f"{sample['box_type']!r} (must be one of {sorted(VALID_BOX_TYPES)})",
                    file=sys.stderr,
                )
                sys.exit(1)
            try:
                distance = float(sample["distance_to_box"])
            except ValueError:
                print(
                    f"Error: sample {sample['sample']!r} has non-numeric "
                    f"distance_to_box {sample['distance_to_box']!r}",
                    file=sys.stderr,
                )
                sys.exit(1)
            if not (0.0 <= distance <= 5.0):
                print(
                    f"Error: sample {sample['sample']!r} distance_to_box {distance} "
                    "must be between 0.0 and 5.0 nm",
                    file=sys.stderr,
                )
                sys.exit(1)

            for key in ("structure", "em_mdp", "nvt_mdp", "npt_mdp", "md_mdp"):
                if not os.path.exists(sample[key]):
                    print(
                        f"Error: sample {sample['sample']!r} {key} file not found: "
                        f"{sample[key]}",
                        file=sys.stderr,
                    )
                    sys.exit(1)

            samples.append(sample)

    if not samples:
        print(f"Error: no samples found in samplesheet {path}", file=sys.stderr)
        sys.exit(1)

    return samples


class GromacsMDWorkflow:
    """Pegasus workflow generator for per-sample GROMACS MD simulations."""

    wf = None
    sc = None
    tc = None
    rc = None
    props = None

    dagfile = None
    wf_dir = None
    shared_scratch_dir = None
    local_storage_dir = None
    wf_name = "gromacs_md_workflow"

    def __init__(self, samples, gmx_cmd="gmx", dagfile="workflow.yml"):
        self.dagfile = dagfile
        self.wf_dir = str(Path(__file__).parent.resolve())
        self.shared_scratch_dir = os.path.join(self.wf_dir, "scratch")
        self.local_storage_dir = os.path.join(self.wf_dir, "output")
        self.samples = samples
        self.gmx_cmd = gmx_cmd

    def write(self):
        if self.sc is not None:
            self.sc.write()
        self.props.write()
        self.rc.write()
        self.tc.write()
        self.wf.write(file=self.dagfile)

    # ------------------------------------------------------------------
    # Plan / run / monitor (thin wrappers over the Pegasus API Workflow
    # object, for interactive use e.g. from a Jupyter notebook)
    # ------------------------------------------------------------------
    def plan_submit(self, exec_site_name="compute"):
        try:
            self.wf.plan(
                dir="submit",
                sites=[exec_site_name],
                output_sites=["local"],
                cleanup="none",
                verbose=1,
                submit=True,
            )
        except PegasusClientError as e:
            print(e)

    def status(self):
        try:
            self.wf.status(long=True)
        except PegasusClientError as e:
            print(e)

    def wait(self):
        try:
            self.wf.wait()
        except PegasusClientError as e:
            print(e)

    def statistics(self):
        try:
            self.wf.statistics()
        except PegasusClientError as e:
            print(e)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    def create_pegasus_properties(self, hosted_site_catalog=None):
        self.props = Properties()
        self.props["pegasus.transfer.threads"] = "16"
        if hosted_site_catalog:
            # Use one of Pegasus' centrally hosted site catalogs instead of
            # a locally generated one. pegasus-plan downloads and caches the
            # named file from the catalog repository at plan time.
            # https://pegasus.isi.edu/documentation/reference-guide/catalogs.html#centrally-hosted-site-catalogs
            self.props["pegasus.catalog.site.repo.file"] = hosted_site_catalog

    # ------------------------------------------------------------------
    # Site Catalog
    #
    # Not used by the CLI below by default — pegasus-plan resolves the site
    # catalog from a centrally hosted one instead (see -s/--hosted-site-catalog
    # and create_pegasus_properties above). Kept for programmatic/notebook use
    # when a self-contained, locally generated HTCondor site catalog is wanted.
    # ------------------------------------------------------------------
    def create_sites_catalog(self, exec_site_name="compute"):
        self.sc = SiteCatalog()

        local = Site("local").add_directories(
            Directory(
                Directory.SHARED_SCRATCH, self.shared_scratch_dir
            ).add_file_servers(
                FileServer("file://" + self.shared_scratch_dir, Operation.ALL)
            ),
            Directory(
                Directory.LOCAL_STORAGE, self.local_storage_dir
            ).add_file_servers(
                FileServer("file://" + self.local_storage_dir, Operation.ALL)
            ),
        )

        exec_site = (
            Site(exec_site_name)
            .add_condor_profile(universe="vanilla")
            .add_pegasus_profile(style="condor")
        )

        self.sc.add_sites(local, exec_site)

    # ------------------------------------------------------------------
    # Transformation Catalog
    # ------------------------------------------------------------------
    def create_transformation_catalog(self, exec_site_name="compute"):
        self.tc = TransformationCatalog()

        container = Container(
            "gromacs_container",
            container_type=Container.SINGULARITY,
            image="file://" + os.path.join(self.wf_dir, "Apptainer", "Gromacs_Container.sif"),
            image_site="local",
        )

        transformations = []
        for tool_name, config in TOOL_CONFIGS.items():
            tx = Transformation(
                tool_name,
                site=exec_site_name,
                pfn=os.path.join(self.wf_dir, f"bin/{tool_name}.py"),
                is_stageable=True,
                container=container,
            ).add_pegasus_profile(
                memory=config["memory"], cores=config.get("cores", 1)
            )
            transformations.append(tx)

        self.tc.add_containers(container)
        self.tc.add_transformations(*transformations)

    # ------------------------------------------------------------------
    # Replica Catalog
    # ------------------------------------------------------------------
    def create_replica_catalog(self):
        self.rc = ReplicaCatalog()

        for s in self.samples:
            sample = s["sample"]
            self.rc.add_replica(
                "local", f"{sample}_structure.pdb", "file://" + os.path.abspath(s["structure"])
            )
            self.rc.add_replica(
                "local", f"{sample}_em.mdp", "file://" + os.path.abspath(s["em_mdp"])
            )
            self.rc.add_replica(
                "local", f"{sample}_nvt.mdp", "file://" + os.path.abspath(s["nvt_mdp"])
            )
            self.rc.add_replica(
                "local", f"{sample}_npt.mdp", "file://" + os.path.abspath(s["npt_mdp"])
            )
            self.rc.add_replica(
                "local", f"{sample}_md.mdp", "file://" + os.path.abspath(s["md_mdp"])
            )

    # ------------------------------------------------------------------
    # Workflow DAG
    # ------------------------------------------------------------------
    def create_workflow(self):
        self.wf = Workflow(self.wf_name, infer_dependencies=True)

        logger.info(f"Creating workflow for {len(self.samples)} sample(s)")
        for s in self.samples:
            self._add_sample_pipeline(s)

    def _add_sample_pipeline(self, s):
        sample = s["sample"]
        gmx = self.gmx_cmd
        logger.info(f"  Adding jobs for sample: {sample}")

        # --- Inputs (Replica Catalog entries) ---
        structure = File(f"{sample}_structure.pdb")
        em_mdp = File(f"{sample}_em.mdp")
        nvt_mdp = File(f"{sample}_nvt.mdp")
        npt_mdp = File(f"{sample}_npt.mdp")
        md_mdp = File(f"{sample}_md.mdp")

        # --- Step 1: pdb_clean_and_check_missing_atoms ---
        checked_pdb = File(f"{sample}_checked.pdb")
        clean_check_job = (
            Job(
                "pdb_clean_and_check_missing_atoms",
                _id=f"pdb_clean_and_check_missing_atoms_{sample}",
                node_label=f"pdb_clean_and_check_missing_atoms_{sample}",
            )
            .add_args("--input", structure, "--output", checked_pdb)
            .add_inputs(structure)
            .add_outputs(checked_pdb, stage_out=False, register_replica=False)
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(clean_check_job)

        # --- Step 2: topology (gmx pdb2gmx) ---
        topo_gro = File(f"{sample}_topology.gro")
        topol_top = File(f"{sample}_topol.top")
        posre_itp = File(f"{sample}_posre.itp")
        topology_job = (
            Job("topology", _id=f"topology_{sample}", node_label=f"topology_{sample}")
            .add_args(
                "--input", checked_pdb,
                "--output-gro", topo_gro,
                "--output-top", topol_top,
                "--output-itp", posre_itp,
                "--force-field", s["force_field"],
                "--gmx-cmd", gmx,
            )
            .add_inputs(checked_pdb)
            .add_outputs(topo_gro, topol_top, posre_itp, stage_out=False, register_replica=False)
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(topology_job)

        # --- Step 3: solvation (editconf + solvate + grompp + genion) ---
        box_solv_ions_gro = File(f"{sample}_box_solv_ions.gro")
        topol_solv_top = File(f"{sample}_topol_solv.top")
        solvation_job = (
            Job("solvation", _id=f"solvation_{sample}", node_label=f"solvation_{sample}")
            .add_args(
                "--input-gro", topo_gro,
                "--input-top", topol_top,
                "--output-gro", box_solv_ions_gro,
                "--output-top", topol_solv_top,
                "--box-type", s["box_type"],
                "--distance-to-box", s["distance_to_box"],
                "--gmx-cmd", gmx,
            )
            .add_inputs(topo_gro, topol_top, posre_itp)
            .add_outputs(
                box_solv_ions_gro, topol_solv_top, stage_out=False, register_replica=False
            )
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(solvation_job)

        # --- Step 4: energy_min (grompp + mdrun) ---
        em_gro = File(f"{sample}_em.gro")
        em_job = (
            Job("energy_min", _id=f"energy_min_{sample}", node_label=f"energy_min_{sample}")
            .add_args(
                "--input-gro", box_solv_ions_gro,
                "--input-top", topol_solv_top,
                "--mdp", em_mdp,
                "--output-gro", em_gro,
                "--gmx-cmd", gmx,
            )
            .add_inputs(box_solv_ions_gro, topol_solv_top, posre_itp, em_mdp)
            .add_outputs(em_gro, stage_out=False, register_replica=False)
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(em_job)

        # --- Step 5: nvt_equilibration (grompp -r + mdrun) ---
        nvt_gro = File(f"{sample}_nvt.gro")
        nvt_job = (
            Job(
                "nvt_equilibration",
                _id=f"nvt_equilibration_{sample}",
                node_label=f"nvt_equilibration_{sample}",
            )
            .add_args(
                "--input-gro", em_gro,
                "--input-top", topol_solv_top,
                "--mdp", nvt_mdp,
                "--output-gro", nvt_gro,
                "--gmx-cmd", gmx,
            )
            .add_inputs(em_gro, topol_solv_top, posre_itp, nvt_mdp)
            .add_outputs(nvt_gro, stage_out=False, register_replica=False)
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(nvt_job)

        # --- Step 6: npt_equilibration (grompp -r + mdrun) ---
        npt_gro = File(f"{sample}_npt.gro")
        npt_job = (
            Job(
                "npt_equilibration",
                _id=f"npt_equilibration_{sample}",
                node_label=f"npt_equilibration_{sample}",
            )
            .add_args(
                "--input-gro", nvt_gro,
                "--input-top", topol_solv_top,
                "--mdp", npt_mdp,
                "--output-gro", npt_gro,
                "--gmx-cmd", gmx,
            )
            .add_inputs(nvt_gro, topol_solv_top, posre_itp, npt_mdp)
            .add_outputs(npt_gro, stage_out=False, register_replica=False)
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(npt_job)

        # --- Step 7: production (grompp + mdrun + report-methods) ---
        md_gro = File(f"{sample}_md.gro")
        md_tpr = File(f"{sample}_md.tpr")
        md_xtc = File(f"{sample}_md.xtc")
        # gmx report-methods requires its -o output to end in .out
        md_report = File(f"{sample}_MD_REPORT.out")
        production_job = (
            Job("production", _id=f"production_{sample}", node_label=f"production_{sample}")
            .add_args(
                "--input-gro", npt_gro,
                "--input-top", topol_solv_top,
                "--mdp", md_mdp,
                "--output-gro", md_gro,
                "--output-tpr", md_tpr,
                "--output-xtc", md_xtc,
                "--output-report", md_report,
                "--gmx-cmd", gmx,
            )
            .add_inputs(npt_gro, topol_solv_top, posre_itp, md_mdp)
            .add_outputs(
                md_gro, md_tpr, md_xtc, md_report,
                stage_out=True, register_replica=False,
            )
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(production_job)

        # --- Step 8: post_processing (trjconv, removes periodicity) ---
        noPBC_xtc = File(f"{sample}_noPBC.xtc")
        post_job = (
            Job(
                "post_processing",
                _id=f"post_processing_{sample}",
                node_label=f"post_processing_{sample}",
            )
            .add_args(
                "--input-tpr", md_tpr,
                "--input-xtc", md_xtc,
                "--output", noPBC_xtc,
                "--gmx-cmd", gmx,
            )
            .add_inputs(md_tpr, md_xtc)
            .add_outputs(noPBC_xtc, stage_out=False, register_replica=False)
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(post_job)

        # --- Step 9: analysis_rmsd (final deliverable) ---
        rmsd_xvg = File(f"{sample}_rmsd.xvg")
        rmsd_job = (
            Job("analysis_rmsd", _id=f"analysis_rmsd_{sample}", node_label=f"analysis_rmsd_{sample}")
            .add_args(
                "--input-gro", md_gro,
                "--input-xtc", noPBC_xtc,
                "--output", rmsd_xvg,
                "--gmx-cmd", gmx,
            )
            .add_inputs(md_gro, noPBC_xtc)
            .add_outputs(rmsd_xvg, stage_out=True, register_replica=False)
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(rmsd_job)

        # --- Step 10: analysis_plot (RMSD .xvg -> PNG; not in the source pipeline) ---
        rmsd_png = File(f"{sample}_rmsd.png")
        plot_job = (
            Job("analysis_plot", _id=f"analysis_plot_{sample}", node_label=f"analysis_plot_{sample}")
            .add_args(
                "--input", rmsd_xvg,
                "--output", rmsd_png,
                "--title", f'"{sample} RMSD vs. Time"',
            )
            .add_inputs(rmsd_xvg)
            .add_outputs(rmsd_png, stage_out=True, register_replica=False)
            .add_pegasus_profiles(label=sample)
        )
        self.wf.add_jobs(plot_job)


def main():
    parser = argparse.ArgumentParser(
        description="Pegasus GROMACS MD Workflow Generator "
        "(converted from nf-core/moleculardynamics)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --samplesheet samplesheet.csv --output workflow.yml
  %(prog)s --samplesheet samplesheet.csv --gmx-cmd gmx_mpi -e compute
  %(prog)s --samplesheet samplesheet.csv -s access-pegasus.yml
""",
    )

    parser.add_argument(
        "--samplesheet",
        metavar="CSV",
        type=str,
        required=True,
        help="CSV with columns: sample,structure,em_mdp,nvt_mdp,npt_mdp,md_mdp,"
        "force_field,box_type,distance_to_box",
    )
    parser.add_argument(
        "--gmx-cmd",
        metavar="STR",
        type=str,
        default="gmx",
        help="GROMACS command to invoke (default: gmx; e.g. gmx_mpi)",
    )
    parser.add_argument(
        "-s",
        "--hosted-site-catalog",
        metavar="FILE",
        type=str,
        default=None,
        help="Name of a Pegasus centrally hosted site catalog to plan against "
        "(e.g. access-pegasus.yml), instead of a locally generated one. Sets "
        "pegasus.catalog.site.repo.file; see "
        "https://pegasus.isi.edu/documentation/reference-guide/catalogs.html"
        "#centrally-hosted-site-catalogs",
    )
    parser.add_argument(
        "-e",
        "--execution-site-name",
        metavar="STR",
        type=str,
        default="compute",
        help="Execution site name (default: compute)",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="STR",
        type=str,
        default="workflow.yml",
        help="Output file (default: workflow.yml)",
    )

    args = parser.parse_args()

    if not os.path.exists(args.samplesheet):
        print(f"Error: samplesheet not found: {args.samplesheet}", file=sys.stderr)
        sys.exit(1)

    samples = parse_samplesheet(args.samplesheet)

    logger.info("=" * 70)
    logger.info("GROMACS MD WORKFLOW GENERATOR")
    logger.info("=" * 70)
    logger.info(f"Samples: {[s['sample'] for s in samples]}")
    logger.info(f"GROMACS command: {args.gmx_cmd}")
    logger.info(f"Execution site: {args.execution_site_name}")
    logger.info(
        f"Hosted site catalog: {args.hosted_site_catalog or '(none — supply your own site catalog)'}"
    )
    logger.info(f"Output file: {args.output}")
    logger.info("=" * 70)

    try:
        workflow = GromacsMDWorkflow(
            samples=samples, gmx_cmd=args.gmx_cmd, dagfile=args.output
        )

        workflow.create_pegasus_properties(hosted_site_catalog=args.hosted_site_catalog)
        workflow.create_transformation_catalog(exec_site_name=args.execution_site_name)
        workflow.create_replica_catalog()
        workflow.create_workflow()
        workflow.write()

        logger.info(f"\nWorkflow written to {args.output}")
        logger.info(
            f"Plan: pegasus-plan --dir submit -s {args.execution_site_name} "
            f"-o local {args.output}"
        )

    except Exception as e:
        logger.error(f"Failed to generate workflow: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
