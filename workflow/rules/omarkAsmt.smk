from os.path import join
import logging
import common_functions as cf


### Configuration Setup ### ----------------------------------------
OMARK_DB = config["omark_database"]


### Logging Setup ### ----------------------------------------
LOG_OMARKASMT = cf.get_logger("OMARKASMT", VERBOSE)


### Rule ### ----------------------------------------
rule omark_assessment:
    """
    Perform Omark assessment to evaluate proteome completeness
    """
    input:
        transcriptPep = join(TRANSEVID_DIR, "transcript.flt.final.pep")
    output:
        omarkOut = directory(join(OMARK_DIR, "omark_output"))
    log:
        join(ASMT_LOG_DIR, "omarkAsmt.log")
    params:
        omarkQuery = join(OMARK_DIR, "query.omamer"),
        omarkDB = OMARK_DB
    threads:
        THREADS
    run:
        LOG_OMARKASMT.info("Running omarkAsmt.smk...")
        startTime = time()

        shell("""
        omamer search \
         --db {params.omarkDB} \
         --query {input.transcriptPep} \
         --out {params.omarkQuery} \
         --nthreads {threads} \
         &> {log}

        omark \
         --file {params.omarkQuery} \
         --database {params.omarkDB} \
         --outputFolder {output.omarkOut} \
         &>> {log}
        """)

        endTime = time()
        elapseTime = endTime - startTime
        LOG_OMARKASMT.info(f"Performed Omark assessment to evaluate proteome completeness in {elapseTime:.2f} seconds")