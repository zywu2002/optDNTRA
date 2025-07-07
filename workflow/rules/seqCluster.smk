from os.path import join
import logging
import optDNTRA_utils as utils
from time import time


### Logging Setup ### ----------------------------------------
LOG_SEQCLUSTER = utils.get_logger("SEQCLUSTER", VERBOSE)


### Rule ### ----------------------------------------
rule cdhit:
    """ 
    Cluster highly similar sequences using CD-HIT-EST
    """
    input:
        transcript=join(PREPROC_DIR, "transcript.fa"),
    output:
        transcriptFlt=join(RMREDUNDANT_DIR, "transcript.flt1.fa"),
    log:
        join(OPTIM_LOG_DIR, "seqCluster.log"),
    threads: THREADS
    run:
        LOG_SEQCLUSTER.info("Running seqCluster.smk...")
        startTime = time()

        shell(
        """
        cd-hit-est \
         -i {input.transcript} \
         -o {output.transcriptFlt} \
         -T {threads} \
         -c 0.95 \
         -M 160000 \
         -n 10 \
         -d 0 \
         -sc 1 \
         -sf 1 \
         &> {log}
        """
        )

        endTime = time()
        elapseTime = endTime - startTime
        LOG_SEQCLUSTER.info(
            f"Clustered highly similar sequences using CD-HIT-EST in {elapseTime:.2f} seconds."
        )
