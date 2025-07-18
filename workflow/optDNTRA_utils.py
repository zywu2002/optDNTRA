#!/usr/bin/env python3

import os
import re
import yaml
import sys
import pandas as pd
import logging
import requests
from tqdm import tqdm
import gzip
import shutil
from coloredlogs import ColoredFormatter
from Bio import SeqIO
import argparse

# Valid file extensions for different file types
fastaExt = ["fa", "fasta"]
fastqExt = ["fq", "fq.gz", "fastq", "fastq.gz"]
hmmExt = ["hmm"]
omarkDbExt = ["h5"]
emapperDbExt = ["db", "db.gz"]

# URLs for reference databases
swiss_prot_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
pfam_hmm_url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
omark_db_url = "https://omabrowser.org/All/LUCA.h5"
emapper_db_url = "http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz"


class OptDNTRAError(Exception):
    """
    Custom exception for optDNTRA workflow errors.
    Raised for any critical error that should halt the pipeline.
    """
    pass


def get_logger(name: str, verbose: bool) -> logging.Logger:
    """
    Create and return a logger with colored output for the given name.

    Args:
        name (str): Logger name.
        verbose (bool): If True, set logging level to DEBUG, else INFO.

    Returns:
        logging.Logger: Configured logger instance.
    """
    logger = logging.getLogger(name)
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    handler = logging.StreamHandler()
    log_formatter = ColoredFormatter("[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    handler.setFormatter(log_formatter)
    logger.addHandler(handler)
    return logger


def load_config(config_file: str, logger: logging.Logger) -> dict:
    """
    Load configuration from a YAML file.

    Args:
        config_file (str): Path to the YAML configuration file.
        logger (logging.Logger): Logger for status and error messages.

    Returns:
        dict: Configuration dictionary loaded from YAML.

    Raises:
        OptDNTRAError: If the file cannot be loaded or parsed.
    """
    logger.info("Loading parameters from the default configuration file...")
    try:
        with open(config_file, "r") as fin:
            config = yaml.load(fin, Loader=yaml.FullLoader)
    except Exception as e:
        logger.error(f"Failed to load config file {config_file}: {e}")
        raise OptDNTRAError(f"Failed to load config file {config_file}: {e}")
    # Replace None values with empty strings for consistency
    config_mod = {k: "" if v is None else v for k, v in config.items()}
    logger.debug("Parameters " + "-" * 110)
    for k, v in config_mod.items():
        logger.debug(f"{k}: {v}")
    logger.debug("-" * 120)
    return config_mod


def load_defaults(args: argparse.Namespace, logger: logging.Logger) -> dict:
    """
    Load the default configuration for the workflow using the provided arguments.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        logger (logging.Logger): Logger for status and error messages.

    Returns:
        dict: Default configuration dictionary.
    """
    config_default = load_config(args.config, logger)
    logger.debug("Loading default configuration for the workflow...")
    return config_default


def check_extensions(file_name: str, ext_list: list[str], logger: logging.Logger) -> None:
    """
    Check if the file exists and has a valid extension.

    Args:
        file_name (str): Path to the file to check.
        ext_list (list[str]): List of valid file extensions.
        logger (logging.Logger): Logger for error messages.

    Raises:
        OptDNTRAError: If the file does not exist or has an invalid extension.
    """
    _, basename = os.path.split(file_name)
    file_extension = basename.split('.', 1)[1] if '.' in basename else ''
    logger.debug(f"Checking the file extension of {file_name}...")
    if not os.path.exists(file_name):
        logger.error(f"Error! Input file {file_name} not found")
        raise OptDNTRAError(f"Input file {file_name} not found")
    if file_extension not in ext_list:
        logger.error(f"Error! {file_name} has an invalid extension")
        raise OptDNTRAError(f"{file_name} has an invalid extension")


def gunzip_file(gzipped_file: str, logger: logging.Logger) -> None:
    """
    Decompress a Gzip file to its original format.

    Args:
        gzipped_file (str): Path to the Gzip file.
        logger (logging.Logger): Logger for error messages.

    Raises:
        OptDNTRAError: If decompression fails.
    """
    logger.debug(f"Decompressing {gzipped_file}...")
    out_file = os.path.splitext(gzipped_file)[0]
    try:
        with gzip.open(gzipped_file, "rb") as fin:
            with open(out_file, "wb") as fout:
                shutil.copyfileobj(fin, fout)
    except Exception as e:
        logger.error(f"Failed to decompress {gzipped_file}: {e}")
        raise OptDNTRAError(f"Failed to decompress {gzipped_file}: {e}")


def download_file(url: str, out_dir: str, logger: logging.Logger) -> None:
    """
    Download a file from a URL and optionally decompress if it is a .gz file.

    Args:
        url (str): URL of the file to download.
        out_dir (str): Output directory for the downloaded file.
        logger (logging.Logger): Logger for status and error messages.

    Raises:
        OptDNTRAError: If the download or decompression fails.
    """
    os.makedirs(out_dir, exist_ok=True)
    file_name = os.path.join(out_dir, url.split("/")[-1])
    logger.info(f"Downloading {file_name} from {url}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))
        block_size = 1024
        with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
            with open(file_name, "wb") as fout:
                for data in response.iter_content(block_size):
                    progress_bar.update(len(data))
                    fout.write(data)
        # Automatically decompress if the file is gzipped
        if file_name.endswith(".gz"):
            gunzip_file(file_name, logger)
    except Exception as e:
        logger.error(f"Failed to download {url}: {e}")
        raise OptDNTRAError(f"Failed to download {url}: {e}")


def load_sampleSheet(file_path: str, config: dict, paired_end: bool = True) -> dict:
    """
    Load a sample sheet from a tab-delimited text file.

    Args:
        file_path (str): Path to the sample sheet file.
        config (dict): Default configuration dictionary.
        paired_end (bool): Whether the data is paired-end (default: True).

    Returns:
        dict: Dictionary mapping sample identifiers to lists of file paths.

    Raises:
        OptDNTRAError: If the sample sheet cannot be read.
    """
    try:
        samples_df = pd.read_csv(file_path, sep="\t", header=None)
    except Exception as e:
        raise OptDNTRAError(f"Failed to read sample sheet {file_path}: {e}")
    if paired_end:
        c2 = samples_df.columns[2]
        c3 = samples_df.columns[3]
        samples_dict = {os.path.basename(row[c2]).split(config["reads"][0])[0]: [row[c2], row[c3]] for _, row in samples_df.iterrows()}
    else:
        c2 = samples_df.columns[2]
        samples_dict = {os.path.basename(row[c2]).split(config["reads"][0])[0]: [row[c2]] for _, row in samples_df.iterrows()}
    return samples_dict


def check_arguments(args: argparse.Namespace, config: dict, logger: logging.Logger) -> None:
    """
    Validate the wrapper arguments and check input file consistency.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        config (dict): Default configuration dictionary.
        logger (logging.Logger): Logger for error messages.

    Raises:
        OptDNTRAError: If argument validation fails or required files are missing.
    """
    logger.info("Checking the wrapper arguments...")
    # Check transcript file
    check_extensions(args.transcript, fastaExt, logger)
    # Check input FASTQ/FASTA logic
    if args.fastq:
        check_extensions(args.fastq, fastqExt, logger)
        if args.left or args.right or args.batch:
            logger.error("Error! invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
            raise OptDNTRAError("Invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
    elif args.left and args.right:
        check_extensions(args.left, fastqExt, logger)
        check_extensions(args.right, fastqExt, logger)
        if args.fastq or args.batch:
            logger.error("Error! invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
            raise OptDNTRAError("Invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
    elif args.batch:
        if args.pairedEnd:
            samples_dict = load_sampleSheet(args.batch, config, paired_end=True)
        else:
            samples_dict = load_sampleSheet(args.batch, config, paired_end=False)
        for _, v in samples_dict.items():
            for f in v:
                check_extensions(f, fastqExt, logger)
        if args.left or args.right or args.fastq:
            logger.error("Error! invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
            raise OptDNTRAError("Invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
    else:
        logger.error("Error! invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
        raise OptDNTRAError("Invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
    # Check reference transcriptome if provided
    if args.reference:
        check_extensions(args.reference, fastaExt, logger)
    # Check and download required databases if missing
    if config["swiss_prot"]:
        check_extensions(config["swiss_prot"], fastaExt, logger)
    else:
        download_file(url=swiss_prot_url, out_dir="data/", logger=logger)
        swiss_prot_path = os.path.splitext(os.path.basename(swiss_prot_url))[0] if swiss_prot_url.endswith(".gz") else os.path.basename(swiss_prot_url)
        config["swiss_prot"] = f"data/{swiss_prot_path}"
    if config["pfam_hmm"]:
        check_extensions(config["pfam_hmm"], hmmExt, logger)
    else:
        download_file(url=pfam_hmm_url, out_dir="data/", logger=logger)
        pfam_hmm_path = os.path.splitext(os.path.basename(pfam_hmm_url))[0] if pfam_hmm_url.endswith(".gz") else os.path.basename(pfam_hmm_url)
        config["pfam_hmm"] = f"data/{pfam_hmm_path}"
    if args.omark:
        if config["omark_database"]:
            check_extensions(config["omark_database"], omarkDbExt, logger)
        else:
            download_file(url=omark_db_url, out_dir="data/", logger=logger)
            omark_db_path = os.path.splitext(os.path.basename(omark_db_url))[0] if omark_db_url.endswith(".gz") else os.path.basename(omark_db_url)
            config["omark_database"] = f"data/{omark_db_path}"
    if args.emapper:
        if config["emapper_database"]:
            check_extensions(config["emapper_database"], emapperDbExt, logger)
        else:
            download_file(url=emapper_db_url, out_dir="data/", logger=logger)
            emapper_db_path = os.path.splitext(os.path.basename(emapper_db_url))[0] if emapper_db_url.endswith(".gz") else os.path.basename(emapper_db_url)
            config["emapper_database"] = f"data/{emapper_db_path}"


def add_sampleConfig(config: dict, logger: logging.Logger) -> dict:
    """
    Add sample information to the configuration dictionary for downstream use.

    Args:
        config (dict): Configuration dictionary.
        logger (logging.Logger): Logger for debug messages.

    Returns:
        dict: Updated configuration dictionary with sample information.
    """
    logger.debug("Adding sample configuration to the config dictionary...")
    if config.get("fastq"):
        config["sample"] = os.path.basename(config["fastq"]).split(".")[0]
    if config.get("left"):
        config["sample"] = os.path.basename(config["left"]).split(config["reads"][0])[0]
    if config.get("batch"):
        config["sample"] = ""
    return config


def write_config(config_file: str, config: dict, logger: logging.Logger) -> None:
    """
    Write the configuration dictionary to a YAML file for Snakemake.

    Args:
        config_file (str): Path to the output YAML file.
        config (dict): Configuration dictionary.
        logger (logging.Logger): Logger for error messages.

    Raises:
        OptDNTRAError: If writing the file fails.
    """
    logger.debug("Writing the configuration dictionary to a YAML file...")
    config = add_sampleConfig(config, logger)
    try:
        with open(config_file, "w") as fin:
            yaml.dump(config, fin, default_flow_style=False, width=1000, sort_keys=False)
    except Exception as e:
        logger.error(f"Failed to write config file {config_file}: {e}")
        raise OptDNTRAError(f"Failed to write config file {config_file}: {e}")


def create_YAML(config_default: dict, args: argparse.Namespace, calling_script: str, logger: logging.Logger) -> None:
    """
    Create a YAML configuration file for Snakemake by merging defaults and command-line arguments.

    Args:
        config_default (dict): Default configuration dictionary.
        args (argparse.Namespace): Parsed command-line arguments.
        calling_script (str): Name of the calling script.
        logger (logging.Logger): Logger for status and debug messages.
    """
    workflow_name = os.path.splitext(os.path.basename(calling_script))[0]
    config_file = os.path.join(f"{workflow_name}.config.yml")
    args_flt_dict = {k: "" if v is None else v for k, v in vars(args).items()}
    del args_flt_dict["snakemakeOptions"]
    config = args_flt_dict
    config.update(config_default)
    write_config(config_file, config, logger)
    logger.info("Creating the YAML configuration file for Snakemake...")
    logger.debug("Parameters " + "-" * 110)
    for k, v in config.items():
        logger.debug(f"{k}: {v}")
    logger.debug("-" * 120)


def count_fastx(fastx: str, logger: logging.Logger, identifier: str = ">") -> int:
    """
    Count the number of sequences in a FASTA or FASTQ file using Biopython.

    Args:
        fastx (str): Path to the FASTA or FASTQ file.
        logger (logging.Logger): Logger for error messages.
        identifier (str): Character identifying the start of a sequence (default is '>').

    Returns:
        int: Number of sequences in the file.

    Raises:
        OptDNTRAError: If counting fails or file is invalid.
    """
    logger.debug(f"Counting the number of sequences in {fastx}...")
    try:
        if identifier == ">":
            seq_count = sum(1 for _ in SeqIO.parse(fastx, "fasta"))
        elif identifier == "@":
            seq_count = sum(1 for _ in SeqIO.parse(fastx, "fastq"))
        else:
            raise ValueError("Unknown identifier for sequence file.")
    except Exception as e:
        logger.error(f"Failed to count sequences in {fastx}: {e}")
        raise OptDNTRAError(f"Failed to count sequences in {fastx}: {e}")
    return seq_count


def update_snakemakeOptions(args: argparse.Namespace, logger: logging.Logger) -> str:
    """
    Update the Snakemake options with additional parameters for execution.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        logger (logging.Logger): Logger for status messages.

    Returns:
        str: A string of snakemake options with '--quiet' and '--cores' parameters added.
    """
    # Add --cores and --quiet (if not verbose) to the snakemake command
    if args.verbose:
        snakemake_cmd = args.snakemakeOptions + [f"--cores {args.threads}"]
    else:
        snakemake_cmd = args.snakemakeOptions + ["--quiet all", f"--cores {args.threads}"]
    snakemake_cmd = " ".join(" ".join(snakemake_cmd).split())
    logger.info(f"Snakemake options: {snakemake_cmd}")
    return snakemake_cmd


def process_path(path: str, new_dir: str) -> str:
    """
    Return a new file path by joining the filename with a new directory.

    Args:
        path (str): Path to the original file.
        new_dir (str): Path to the new directory.

    Returns:
        str: New file path in the new directory.
    """
    filename = os.path.basename(path)
    return os.path.join(new_dir, filename)
