from os.path import join
from os import makedirs
import logging
import optDNTRA_utils as utils
from time import time


### Configuration Setup ### ----------------------------------------
FASTQ = config["fastq"]
FASTQ_LEFT = config["left"]
FASTQ_RIGHT = config["right"]


### Logging Setup ### ----------------------------------------
LOG_QUALCTRL = utils.get_logger("QUALITYCTRL", VERBOSE)


### Rule ### ----------------------------------------
if BATCH:
    if PAIREDEND:

        rule fastqc:
            """
            Perform quality check using FastQC
            """
            input:
                fastqCleanLeft=join(PREPROC_DIR, "{sample}" + READS[0] + FASTQ_EXT),
                fastqCleanRight=join(PREPROC_DIR, "{sample}" + READS[1] + FASTQ_EXT),
            output:
                fastqcOutDir=directory(join(FASTQC_DIR, "{sample}")),
            log:
                join(PREPROC_LOG_DIR, "qualityCtrl.{sample}.log"),
            threads: THREADS
            run:
                LOG_QUALCTRL.info("Running qualityCtrl.smk...")
                startTime = time()

                makedirs(output.fastqcOutDir, exist_ok=True)

                shell(
                """
                fastqc \
                --threads {threads} \
                --outdir {output.fastqcOutDir} \
                {input.fastqCleanLeft} \
                {input.fastqCleanRight} \
                &> {log}
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_QUALCTRL.info(
                    f"Performed quality check of {wildcards.sample} using FastQC in {elapseTime:.2f} seconds."
                )

    else:

        rule fastqc:
            """
            Perform quality check using FastQC
            """
            input:
                fastqClean=join(PREPROC_DIR, "{sample}" + FASTQ_EXT),
            output:
                fastqcOutDir=directory(join(FASTQC_DIR, "{sample}")),
            log:
                join(PREPROC_LOG_DIR, "qualityCtrl.{sample}.log"),
            threads: THREADS
            run:
                LOG_QUALCTRL.info("Running qualityCtrl.smk...")
                startTime = time()

                makedirs(output.fastqcOutDir, exist_ok=True)

                shell(
                """
                fastqc \
                --threads {threads} \
                --outdir {output.fastqcOutDir} \
                {input.fastqClean} \
                &> {log}
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_QUALCTRL.info(
                    f"Performed quality check of {wildcards.sample} using FastQC in {elapseTime:.2f} seconds."
                )

else:
    if PAIREDEND:

        rule fastqc:
            """
            Perform quality check using FastQC
            """
            input:
                fastqCleanLeft=join(PREPROC_DIR, SAMPLE + READS[0] + FASTQ_EXT),
                fastqCleanRight=join(PREPROC_DIR, SAMPLE + READS[1] + FASTQ_EXT),
            output:
                fastqcOutDir=directory(FASTQC_DIR),
            log:
                join(PREPROC_LOG_DIR, "qualityCtrl.log"),
            threads: THREADS
            run:
                LOG_QUALCTRL.info("Running qualityCtrl.smk...")
                startTime = time()

                makedirs(output.fastqcOutDir, exist_ok=True)

                shell(
                """
                fastqc \
                --threads {threads} \
                --outdir {output.fastqcOutDir} \
                {input.fastqCleanLeft} \
                {input.fastqCleanRight} \
                &> {log}
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_QUALCTRL.info(
                    f"Performed quality check using FastQC in {elapseTime:.2f} seconds."
                )

    else:

        rule fastqc:
            """
            Perform quality check using FastQC
            """
            input:
                fastqClean=join(PREPROC_DIR, SAMPLE + FASTQ_EXT),
            output:
                fastqcOutDir=directory(FASTQC_DIR),
            log:
                join(PREPROC_LOG_DIR, "qualityCtrl.log"),
            threads: THREADS
            run:
                LOG_QUALCTRL.info("Running qualityCtrl.smk...")
                startTime = time()

                makedirs(output.fastqcOutDir, exist_ok=True)

                shell(
                """
                fastqc \
                --threads {threads} \
                --outdir {output.fastqcOutDir} \
                {input.fastqClean} \
                &> {log}
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_QUALCTRL.info(
                    f"Performed quality check using FastQC in {elapseTime:.2f} seconds."
                )
