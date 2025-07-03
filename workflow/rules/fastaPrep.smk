from os.path import join
import logging
import common_functions as cf
from time import time


### Configuration Setup ### ----------------------------------------
TRANSCRIPT = config["transcript"]


### Logging Setup ### ----------------------------------------------
LOG_FASTAPREP = cf.get_logger("FASTAPREP", VERBOSE)


### Rule ### -------------------------------------------------------
rule prepare_fasta:
    """
    Preprocess input FASTA file
    """
    input:
        transcript=TRANSCRIPT,
    output:
        transcriptProc=join(PREPROC_DIR, "transcript.fa"),
    log:
        join(PREPROC_LOG_DIR, "fastaPrep.log"),
    threads: THREADS
    run:
        LOG_FASTAPREP.info("Running fastaPrep.smk...")
        startTime = time()

        shell(
        """
        sed 's/>\\([^ ]*\\) .*/>\\1/' {input.transcript} > {output.transcriptProc} 2> {log}
        """
        )

        endTime = time()
        elapseTime = endTime - startTime
        LOG_FASTAPREP.info(
            f"Preprocessed input FASTA file in {elapseTime:.2f} seconds."
        )
