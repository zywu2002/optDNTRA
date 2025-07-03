from os.path import join
from os import makedirs
import logging
import common_functions as cf
from time import time


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
        dataDir=directory(join(ANNO_DIR, "database")),
    threads: THREADS
    run:
        LOG_EMAPPER.info("Running emapper.smk...")
        startTime = time()

        makedirs(params.dataDir, exist_ok=True)

        shell(
        """
        download_eggnog_data.py -y --data_dir {params.dataDir}

        emapper.py \
         -m diamond \
         --itype proteins \
         -i {input.transcriptPep} \
         --data_dir {params.dataDir} \
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
