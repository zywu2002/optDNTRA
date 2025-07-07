#!/usr/bin/env python3

import argparse
import sys
import workflow.optDNTRA_utils as utils
import subprocess
from time import time
import os


def parse_args() -> argparse.ArgumentParser:
    """
    Parse command-line arguments for the optDNTRA pipeline.

    Returns:
        argparse.ArgumentParser: Configured argument parser instance.
    """
    parser = argparse.ArgumentParser(
        prog="optDNTRA.py",
        description='''
optDNTRA: Optimization of De Novo Transcriptome RNA-seq Assembly

For RNA-seq input data:
    If Paired-end:
        -1 <string>, --left <string>         left reads
        -2 <string>, --right <string>        right reads
    Or Single-end:
        -f <string>, --fastq <string>        single-end reads
    Or
        -s <string>, --sampleSheet <string>  tab-delimited text file indicating biological replicate relationships
        Example:
            cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
            cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
            cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
            cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
        For single-end, remove the 4th column in the text file.
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Thank you for using optDNTRA (Optimization of De Novo Transcriptome RNA-seq Assembly)!'
    )

    # Workflow and input arguments
    parser.add_argument("-c", "--config",
                        type=str,
                        required=True,
                        help="Path to the workflow configuration file (defaults.yml)",
                        metavar="defaults.yml",
                        dest="config")

    parser.add_argument("-t", "--transcript",
                        type=str,
                        required=True,
                        help="Transcript FASTA file",
                        metavar="transcripts.fasta",
                        dest="transcript")

    parser.add_argument("-f", "--fastq",
                        type=str,
                        required=False,
                        help="Single-end reads (FASTQ)",
                        metavar="reads.fq",
                        dest="fastq")

    parser.add_argument("-1", "--left",
                        type=str,
                        required=False,
                        help="Left reads for paired-end data (FASTQ)",
                        metavar="reads_1.fq",
                        dest="left")

    parser.add_argument("-2", "--right",
                        type=str,
                        required=False,
                        help="Right reads for paired-end data (FASTQ)",
                        metavar="reads_2.fq",
                        dest="right")

    parser.add_argument("-s", "--sampleSheet",
                        type=str,
                        required=False,
                        help="Tab-delimited file indicating biological replicate relationships",
                        metavar="samples.tab",
                        dest="batch")

    parser.add_argument("-o", "--outDir",
                        type=str,
                        default="optDNTRA_out",
                        required=False,
                        help="Output directory (default: optDNTRA_out)",
                        metavar="optDNTRA_out",
                        dest="outDir")

    parser.add_argument("-r", "--reference",
                        type=str,
                        required=False,
                        help="Reference transcriptome (FASTA)",
                        metavar="reference.fasta",
                        dest="reference")

    parser.add_argument("-se", "--singleEnd",
                        action="store_false",
                        help="Specify if the input data is single-end (not paired-end)",
                        dest="pairedEnd")

    parser.add_argument("-ss", "--ss-lib-type",
                        type=str,
                        required=False,
                        choices=["F", "R", "RF", "FR"],
                        help="Strand-specific library type: single ('F' or 'R'), paired ('RF' or 'FR')",
                        dest="ssLibType")

    parser.add_argument("--trim",
                        action="store_true",
                        help="Enable trimming for input data",
                        dest="trim")

    parser.add_argument("--qc",
                        action="store_true",
                        help="Enable quality control for input data",
                        dest="qc")

    parser.add_argument("--buscoAsmt",
                        action="store_true",
                        help="Enable BUSCO assessment",
                        dest="busco")

    parser.add_argument("--omarkAsmt",
                        action="store_true",
                        help="Enable OMArk assessment",
                        dest="omark")

    parser.add_argument("--emapperAnno",
                        action="store_true",
                        help="Enable EggNOG-mapper for functional annotation",
                        dest="emapper")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Print detailed reports (verbose mode)",
                        dest="verbose")

    parser.add_argument("-p", "--threads",
                        type=int,
                        default=1,
                        required=False,
                        help="Number of threads to use (default: 1)",
                        dest="threads")

    parser.add_argument("--snakemakeOptions",
                        type=str,
                        required=False,
                        help="Additional options to pass directly to Snakemake, e.g. --snakemakeOptions='--dryrun'",
                        dest="snakemakeOptions")

    return parser


def main() -> None:
    """
    Main entry point for the optDNTRA workflow.
    Handles argument parsing, configuration loading, pipeline validation, and Snakemake execution.
    """
    main_start_time = time()
    parser = parse_args()
    args = parser.parse_args()
    LOG_MAIN = utils.get_logger("MAIN", getattr(args, "verbose", False))
    try:
        # Always resolve the absolute path to the workflow/Snakefile relative to this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        snakefile_path = os.path.join(script_dir, "workflow", "Snakefile")
        # Always prepend the correct snakemake --snakefile path to the options
        # If user provided additional snakemake options, append them; otherwise, only use the default
        args.snakemakeOptions = [f"snakemake --snakefile {snakefile_path}"] + [args.snakemakeOptions] if args.snakemakeOptions else [f"snakemake --snakefile {snakefile_path}"]

        # Load workflow configuration (defaults.yml and command-line args)
        config_default = utils.load_defaults(args, LOG_MAIN)
        # Validate all arguments and configuration
        utils.check_arguments(args, config_default, LOG_MAIN)
        # Generate the YAML configuration file for Snakemake
        utils.create_YAML(config_default, args, __file__, LOG_MAIN)
        # Build the full Snakemake command string
        snakemake_cmd = utils.update_snakemakeOptions(args, LOG_MAIN)
        # Run the Snakemake workflow
        LOG_MAIN.info("Running Snakemake...")
        subprocess.run(snakemake_cmd, shell=True, check=True)
        # Log total runtime
        main_end_time = time()
        elapse_time = main_end_time - main_start_time
        LOG_MAIN.info(f"Total runtime of the project: {elapse_time:.2f} seconds")
    except utils.OptDNTRAError as e:
        LOG_MAIN.error(f"optDNTRA workflow error: {e}")
        sys.exit(1)
    except Exception as e:
        LOG_MAIN.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
