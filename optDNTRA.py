#!/usr/bin/env python3

import argparse
import sys
import workflow.common_functions as cf
import subprocess
from time import time


def parse_args():
    """ Parse arguments from the command line

    Returns:
        parser (argparse.ArgumentParser): Parsed arguments
    """
    # Parser arguments
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description=''' 
optDNTRA: optimization of De Novo Transcriptome Rna-seq Assembly

For RNA-seq input data:
    If Paired-end:
        --left <string>         left reads
        --right <string>        right reads
    or Single-end:
        --fastq, -f <string>    single-end reads
    or
        --sampleSheet <string>  tab-delimited text file indicating biological replicate relationships
                                e.g.
                                cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
                                cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
                                cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
                                cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
                                if single-end, remove the 4th column in the text file.
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Thank you for using optDNTRA!'
    )
    
    # Define workflow arguments
    parser.add_argument("-t", "--transcript",
                        type=str,
                        required=True,
                        help="transcript fasta",
                        metavar="transcripts.fasta",
                        dest="transcript")
    
    parser.add_argument("-f", "--fastq",
                        type=str,
                        required=False,
                        help="single-end reads",
                        metavar="reads.fq",
                        dest="fastq")

    parser.add_argument("--left",
                        type=str,
                        required=False,
                        help="left reads",
                        metavar="reads_1.fq",
                        dest="left")
    
    parser.add_argument("--right",
                        type=str,
                        required=False,
                        help="right reads",
                        metavar="reads_2.fq",
                        dest="right")
    
    parser.add_argument("--sampleSheet",
                        type=str,
                        required=False,
                        help="tab-delimited text file indicating biological replicate relationships",
                        metavar="samples.tab",
                        dest="batch")
    
    parser.add_argument("--reference",
                        type=str,
                        required=False,
                        help="reference transcripts",
                        metavar="reference.fasta",
                        dest="reference")
    
    parser.add_argument("--singleEnd",
                        action="store_false",
                        help="specify if the input data is single-end, not paired-end",
                        dest="pairedEnd")
    
    parser.add_argument("--ss-lib-type",
                        type=str,
                        required=False,
                        choices=["F", "R", "RF", "FR"],
                        help="strand-specific library type: single('F' or 'R'), paired('RF' or 'FR')",
                        dest="ssLibType")
    
    parser.add_argument("--trim",
                        action="store_true",
                        help="enable trimming for input data",
                        dest="trim")
    
    parser.add_argument("--qc",
                        action="store_true",
                        help="enable quality control for input data",
                        dest="qc")

    parser.add_argument("--buscoAsmt",
                        action="store_true",
                        help="enable BUSCO assessment",
                        dest="busco")
    
    parser.add_argument("--omarkAsmt",
                        action="store_true",
                        help="enable OMArk assessment",
                        dest="omark")

    parser.add_argument("--emapper",
                        action="store_true",
                        help="enable EggNOG-mapper for functional annotation",
                        dest="emapper")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="print detailed reports",
                        dest="verbose")
    
    parser.add_argument("--threads",
                        type=int,
                        default=1,
                        required=False,
                        help="number of threads to use",
                        dest="threads")
    
    parser.add_argument("--snakemakeOptions",
                        action="append",
                        default=["snakemake --snakefile workflow/Snakefile"],
                        help="Snakemake options to be passed directly to snakemake, e.g. use --snakemakeOptions='--dryrun'",
                        dest="snakemakeOptions")

    return parser


def main():
    mainStartTime = time()
    
    # Parse command-line arguments
    parser = parse_args()
    args = parser.parse_args()
    
    # Configure the logger
    LOG_MAIN = cf.get_logger("MAIN", args.verbose)
    
    # Load default configurations
    baseDir, workflowDir, configDefault = cf.load_defaults(__file__, LOG_MAIN)

    # Check and validate the wrapper arguments
    cf.check_arguments(args, configDefault, LOG_MAIN)
    
    # Generate the YAML configuration file for Snakemake
    cf.create_YAML(configDefault, args, __file__, LOG_MAIN)
    
    # Update Snakemake options
    snakemakeCMD = cf.update_snakemakeOptions(args, LOG_MAIN)
    
    # Execute Snakemake
    LOG_MAIN.info("Running Snakemake...")
    subprocess.run(snakemakeCMD, shell=True)
    
    # Logging out runtime
    mainEndTime = time()
    elapseTime = mainEndTime - mainStartTime
    LOG_MAIN.info(f"Total runtime of the project: {elapseTime:.2f} seconds")


if __name__ == "__main__":
    main()
