from os.path import join, dirname
from os import makedirs
import logging
import common_functions as cf
from time import time


### Configuration Setup ### ----------------------------------------
EMAPPER_DB = dirname(config["emapper_database"])


### Logging Setup ### ----------------------------------------
LOG_EMAPPER = cf.get_logger("EMAPPER", VERBOSE)


### Rule ### ----------------------------------------
rule emapper:
    """
    Perform EggNOG-mapper for functional annotation
    """
    input:
        transcriptPep=join(TRANSEVID_DIR, "transcript.flt.final.pep"),
    output:
        emapperOut=directory(ANNO_DIR),
    log:
        join(ANNO_LOG_DIR, "emapper.log"),
    params:
        emapperDB=EMAPPER_DB,
    threads: THREADS
    run:
        LOG_EMAPPER.info("Running emapper.smk...")
        startTime = time()

        makedirs(output.emapperOut, exist_ok=True)

        shell(
        """
        emapper.py \
         -m diamond \
         --itype proteins \
         -i {input.transcriptPep} \
         --data_dir {params.emapperDB} \
         --output transAsm \
         --output_dir {output.emapperOut} \
         --cpu {threads} \
         &> {log}
         """
        )

        endTime = time()
        elapseTime = endTime - startTime
        LOG_EMAPPER.info(
            f"Performed EggNOG-mapper for functional annotation in {elapseTime:.2f} seconds."
        )
