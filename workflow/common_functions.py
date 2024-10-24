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


# Valid extensions
fastaExt = ["fa", "fasta"]
fastqExt = ["fq", "fq.gz", "fastq", "fastq.gz"]
hmmExt = ["hmm"]
omarkDbExt = ["h5"]

# Database URL
swiss_prot_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
pfam_hmm_url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
omark_db_url = "https://omabrowser.org/All/LUCA.h5"

# Functions
def get_logger(name, verbose):
    """ Return a logger with a default ColoredFormatter

    Args:
        name (str): logger name
        verbose (bool): If True, set logging level to DEBUG

    Returns:
        logger (logging.Logger): Configured logger instance
    """
    logger = logging.getLogger(name)
    if logger.hasHandlers:
        logger.handlers.clear()
    logger.setLevel(logging.DEBUG) if verbose else logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    logFormatter = ColoredFormatter("[%(asctime)s] %(levelname)s %(name)s: %(message)s")
    handler.setFormatter(logFormatter)
    logger.addHandler(handler)
    return logger


def load_config(configFile, logger):
    """ Load configuration from a YAML file

    Args:
        configFile (str): Path to the YAML configuration file
        logger (logging.Logger): Configured logger instance

    Returns:
        config (dict): Configuration dictionary
    """
    logger.info("Loading parameters from the default configuration file...")
    
    with open(configFile, "r") as FIN:
        config = yaml.load(FIN, Loader=yaml.FullLoader)
    
    configMod = {k: "" if v is None else v for k, v in config.items()}
    
    logger.debug("**Parameters** " + "-" * 100)
    for k, v in configMod.items():
        logger.debug(f"{k}: {v}")
    logger.debug("-" * 120)
    return configMod


def load_defaults(fileName, logger):
    """ Load default configuration for the workflow

    Args:
        fileName (str): Path to the wrapper file
        logger (logging.Logger): Configured logger instance

    Returns:
        baseDir (str): Base directory
        configDir (str): Configuration directory
        workflowDir (str): Workflow directory
        configDefault (dict): Default configuration dictionary
    """
    baseDir = os.path.dirname(fileName)
    configDir = os.path.join(baseDir, "config")
    workflowDir = os.path.join(baseDir, "workflow")

    configFile = os.path.join(configDir, "defaults.yaml")
    configDefault = load_config(configFile, logger)

    logger.debug("Loading default configuration for the workflow...")
    return baseDir, configDir, workflowDir, configDefault


def check_extensions(fileName, extLst, logger):
    """ Check the file extension against a list of valid extensions

    Args:
        fileName (str): Path to the file
        extLst (list): List of valid extensions
    """
    _, basename = os.path.split(fileName)
    fileExtension = basename.split('.', 1)[1]
    
    logger.debug(f"Checking the file extension of {fileName}...")
    
    if not os.path.exists(fileName):
        logger.error(f"Error! Input file {fileName} not found")
        sys.exit(1)
    
    if fileExtension not in extLst:
        logger.error(f"Error! {fileName} has an invalid extension")
        sys.exit(1)


def gunzip_file(gzippedFile, logger):
    """ Decompress a Gzip file

    Args:
        gzippedFile (str): Path to the Gzip file
    """
    logger.debug(f"Decompressing {gzippedFile}...")
    
    outFile = os.path.splitext(gzippedFile)[0]
    with gzip.open(gzippedFile, "rb") as FIN:
        with open(outFile, "wb") as FOUT:
            shutil.copyfileobj(FIN, FOUT)


def download_file(url, outDir, logger):
    """ Download a file from a URL

    Args:
        url (str): URL of a file to download
        outDir (str): Output directory for the downloaded file
    """
    logger.info(f"Downloading {fileName} from {url}...")
    
    os.makedirs(outDir, exist_ok=True)
    fileName = os.path.join(outDir, url.split("/")[-1])
    
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get("content-length", 0))
    block_size = 1024
    with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
        with open(fileName, "wb") as FOUT:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                FOUT.write(data)
        
    if fileName.endswith(".gz"):
        gunzip_file(fileName, logger)


def load_sampleSheet(filePath, config, pairedEnd=True):
    """ Load samples sheet from a text file

    Args:
        filePath (str): Path to the text file containing sample information
        config (dict): Default configuration dictionary

    Returns:
        samplesDict (dict): A dictionary where keys are sample identifiers and values are lists containing sample paths
    """
    
    samplesDf = pd.read_csv(filePath, sep="\t", header=None)
    if pairedEnd:
        c2 = samplesDf.columns[2]
        c3 = samplesDf.columns[3]
        samplesDict = {os.path.basename(row[c2]).split(config["reads"][0])[0]: [row[c2], row[c3]] for _, row in samplesDf.iterrows()}
    else:
        c2 = samplesDf.columns[2]
        samplesDict = {os.path.basename(row[c2]).split(config["reads"][0])[0]: [row[c2]] for _, row in samplesDf.iterrows()}
    return samplesDict


def check_arguments(args, config, logger):
    """ Validate the wrapper arguments

    Args:
        args (argparse.Namespace): Wrapper arguments
        config (dict): Default configuration dictionary
    """
    logger.info("Checking the wrapper arguments...")
    
    check_extensions(args.transcript, fastaExt, logger)
    
    if args.fastq:
        check_extensions(args.fastq, fastqExt, logger)
        if args.left or args.right or args.batch:
            logger.error("Error! invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
            sys.exit(1)
            
    elif args.left and args.right:
        check_extensions(args.left, fastqExt, logger)
        check_extensions(args.right, fastqExt, logger)
        if args.fastq or args.batch:
            logger.error("Error! invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
            sys.exit(1)
            
    elif args.batch:
        if args.pairedEnd:
            samplesDict = load_sampleSheet(args.batch, config, pairedEnd=True)
        else:
            samplesDict = load_sampleSheet(args.batch, config, pairedEnd=False)
        for k, v in samplesDict.items():
            for f in v:
                check_extensions(f, fastqExt, logger)
        if args.left or args.right or args.fastq:
            logger.error("Error! invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
            sys.exit(1)
            
    else:
        logger.error("Error! invalid fastq input! please use ['-f', '--fastq'] or ['--left', '--right'] or ['-b', '--batch']")
        sys.exit(1)
    
    if args.reference:
        check_extensions(args.reference, fastaExt, logger)
        
    if config["swiss_prot"]:
        check_extensions(config["swiss_prot"], fastaExt, logger)
    else:
        download_file(url=swiss_prot_url, outDir="data/", logger=logger)
        
    if config["pfam_hmm"]:
        check_extensions(config["pfam_hmm"], hmmExt, logger)
    else:
        download_file(url=pfam_hmm_url, outDir="data/", logger=logger)
        
    if config["omark_database"]:
        check_extensions(config["omark_database"], omarkDbExt, logger)
    else:
        download_file(url=omark_db_url, outDir="data/", logger=logger)


def add_sampleConfig(config, logger):
    """ Add sample configuration to the config dictionary

    Args:
        config (dict): Configuration dictionary

    Returns:
        config (dict): Updated configuration dictionary
    """
    logger.debug("Adding sample configuration to the config dictionary...")
    
    if config["fastq"]:
        config["sample"] = os.path.basename(config["fastq"]).split(".")[0]
    if config["left"] and config["left"]:
        config["sample"] = os.path.basename(config["left"]).split(config["reads"][0])[0]
    if config["batch"]:
        config["sample"] = ""
    return config


def write_config(configFile, config, logger):
    """ Write the configuration dictionary to a YAML file

    Args:
        configFile (str): Path to the output configuration file
        config (dict): Configuration dictionary
    """
    logger.debug("Writing the configuration dictionary to a YAML file...")
    
    config = add_sampleConfig(config, logger)
    with open(configFile, "w") as FIN:
        yaml.dump(config, FIN, default_flow_style=False, width=1000, sort_keys=False)


def create_YAML(configDir, configDefault, args, callingScript, logger):
    """ Create a YAML configuration file for Snakemake

    Args:
        configDir (str): Directory of the configuration file
        configDefault (dict): Default configuration dictionary
        args (argparse.Namespace): Parsed arguments
        callingScript (str): Name of the calling script
    """
    workflowName = os.path.splitext(os.path.basename(callingScript))[0]
    configFile = os.path.join(configDir, f"{workflowName}.config.yaml")

    argsFltDict = {k: "" if v is None else v for k, v in vars(args).items()}
    del(argsFltDict["snakemakeOptions"])
    config = argsFltDict
    config.update(configDefault)
    write_config(configFile, config, logger)
    
    logger.info("Creating the YAML configuration file for Snakemake...")
    logger.debug("**Parameters** " + "-" * 100)
    for k, v in config.items():
        logger.debug(f"{k}: {v}")
    logger.debug("-" * 120)
    
    
def count_fastx(fastx, logger, identifier=">"):
    """ Count the number of sequences in a FASTA or FASTQ file

    Args:
        fastx (str): Path to the FASTA or FASTQ file
        identifier (str): Character identifying the start of a sequence (default is '>')

    Returns:
        seqCount (int): Number of sequences in the file
    """
    logger.debug(f"Counting the number of sequences in a {fastx}...")
    
    if identifier == ">":
        seqCount = sum(1 for _ in SeqIO.parse(fastx, "fasta"))
    if identifier == "@":
        seqCount = sum(1 for _ in SeqIO.parse(fastx, "fastq"))
    return seqCount
    
    
def update_snakemakeOptions(args, logger):
    """ Update the Snakemake options with additional parameters
    
    Args:
        args (argparse.Namespace): Wrapper arguments

    Returns:
        snakemakeCMD (str): A string of snakemake options with '--quiet' and '--cores' parameters added
    """
    if args.verbose:
        snakemakeCMD = args.snakemakeOptions + [f"--cores {args.threads}"]
    else:
        snakemakeCMD = args.snakemakeOptions + ["--quiet all", f"--cores {args.threads}"]
    snakemakeCMD = " ".join(" ".join(snakemakeCMD).split())
    
    logger.info(f"Snakemake options: {snakemakeCMD}")
    return snakemakeCMD
    
    
def process_path(path, newDir):
    """ Process file path with new directory
    
    Args:
        path (str): Path to the original file
        newDir (str): Path to the new directory
    """
    filename = os.path.basename(path)
    return os.path.join(newDir, filename)
