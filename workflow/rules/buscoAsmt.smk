from os.path import join
import logging
import optDNTRA_utils as utils
from time import time


### Configuration Setup ### ----------------------------------------
TRANSCRIPTPREFLT = config["transcript"]


### Logging Setup ### ----------------------------------------
LOG_BUSCOASMT = utils.get_logger("BUSCOASMT", VERBOSE)


### Rule ### ----------------------------------------
checkpoint busco_assessment:
    """
    Perform BUSCO assessment to evaluate transcriptome completeness
    """
    input:
        transcriptPreFlt=TRANSCRIPTPREFLT,
        transcriptPostFlt=join(TRANSEVID_DIR, "transcript.flt.final.fa"),
    output:
        buscoPreFlt=directory(
            join(BUSCO_DIR, BUSCO_LINEAGE, "busco_" + BUSCO_LINEAGE + "_preFlt")
        ),
        buscoPostFlt=directory(
            join(BUSCO_DIR, BUSCO_LINEAGE, "busco_" + BUSCO_LINEAGE + "_postFlt")
        ),
    log:
        join(ASMT_LOG_DIR, "buscoAsmt-busco_assessment.log"),
    params:
        lineage=BUSCO_LINEAGE,
        buscoPreFlt="busco_" + BUSCO_LINEAGE + "_preFlt",
        buscoPostFlt="busco_" + BUSCO_LINEAGE + "_postFlt",
        downloadPath=directory(join(BUSCO_DIR, BUSCO_LINEAGE))
    threads: THREADS
    run:
        LOG_BUSCOASMT.info("Running buscoAsmt.smk...")
        startTime = time()

        shell(
        """
        busco \
         --in {input.transcriptPreFlt} \
         --out {params.buscoPreFlt} \
         --mode transcriptome \
         --lineage_dataset {params.lineage} \
         --cpu {threads} \
         --download_path {params.downloadPath} \
         --force \
         --out_path {params.downloadPath} \
         --quiet \
         &> {log}

        busco \
         --in {input.transcriptPostFlt} \
         --out {params.buscoPostFlt} \
         --mode transcriptome \
         --lineage_dataset {params.lineage} \
         --cpu {threads} \
         --download_path {params.downloadPath} \
         --force \
         --out_path {params.downloadPath} \
         --quiet \
         &>> {log}
        """
        )

        endTime = time()
        elapseTime = endTime - startTime
        LOG_BUSCOASMT.info(
            f"Performed BUSCO assessment for transcriptome completeness in {elapseTime:.2f} seconds."
        )


def get_busco_summary(wildcards):
    buscoPreFlt = checkpoints.busco_assessment.get(**wildcards).output.buscoPreFlt
    buscoPostFlt = checkpoints.busco_assessment.get(**wildcards).output.buscoPostFlt
    checkpointOut = [
        join(
            buscoPreFlt,
            f"short_summary.specific.{BUSCO_LINEAGE}.busco_{BUSCO_LINEAGE}_preFlt.txt",
        ),
        join(
            buscoPostFlt,
            f"short_summary.specific.{BUSCO_LINEAGE}.busco_{BUSCO_LINEAGE}_postFlt.txt",
        ),
    ]
    return checkpointOut


rule generate_plot:
    """
    Generate BUSCO assessment plot
    """
    input:
        buscoSummary=get_busco_summary,
    output:
        buscoFig=join(BUSCO_DIR, BUSCO_LINEAGE, "busco_figure.png"),
    log:
        join(ASMT_LOG_DIR, "buscoAsmt-generate_plot.log"),
    params:
        buscoDir=directory(join(BUSCO_DIR, BUSCO_LINEAGE)),
        plotScript="generate_plot.py",
    threads: THREADS
    run:
        shell(
        """
        cp {input.buscoSummary[0]} {params.buscoDir} &> {log}
        cp {input.buscoSummary[1]} {params.buscoDir} &>> {log}
        {params.plotScript} \
         --working_directory {params.buscoDir} \
         --quiet \
         &>> {log}
        touch {output.buscoFig}
        """
        )
