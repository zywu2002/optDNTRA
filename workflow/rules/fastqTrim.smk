from os.path import join, splitext, basename, dirname, exists
from os import makedirs
import logging
import optDNTRA_utils as utils
from time import time
import shutil
import pandas as pd


### Configuration Setup ### ----------------------------------------
RAWDATA_PATH = "data"


### Logging Setup ### ----------------------------------------------
LOG_FQTRIM = utils.get_logger("FASTQTRIM", VERBOSE)


### Rule ### -------------------------------------------------------
if TRIM:
    if BATCH:
        if PAIREDEND:

            rule fastp_cleanFQ:
                """
                Preprocess input FASTQ file
                """
                input:
                    fastqLeft=lambda wildcards: SAMPLE_DICT[wildcards.sample][0],
                    fastqRight=lambda wildcards: SAMPLE_DICT[wildcards.sample][1],
                output:
                    fastqCleanLeft=join(PREPROC_DIR, "{sample}" + READS[0] + FASTQ_EXT),
                    fastqCleanRight=join(PREPROC_DIR, "{sample}" + READS[1] + FASTQ_EXT),
                    fastpHTML=join(PREPROC_DIR, "fastp.{sample}.html"),
                    fastpJSON=join(PREPROC_DIR, "fastp.{sample}.json"),
                log:
                    join(PREPROC_LOG_DIR, "fastqTrim.{sample}.log"),
                threads: THREADS
                run:
                    LOG_FQTRIM.info("Running fastqTrim.smk...")
                    startTime = time()

                    shell(
                    """
                    fastp \
                     --thread {threads} \
                     --in1 {input.fastqLeft} --in2 {input.fastqRight} \
                     --out1 {output.fastqCleanLeft} --out2 {output.fastqCleanRight} \
                     --html {output.fastpHTML} \
                     --json {output.fastpJSON} \
                     &> {log}
                    """
                    )

                    endTime = time()
                    elapseTime = endTime - startTime
                    LOG_FQTRIM.info(
                        f"Preprocessed input FASTQ file of {wildcards.sample} in {elapseTime:.2f} seconds"
                    )

        else:

            rule fastp_cleanFQ:
                """
                Preprocess input FASTQ file
                """
                input:
                    fastq=lambda wildcards: SAMPLE_DICT[wildcards.sample][0],
                output:
                    fastqClean=join(PREPROC_DIR, "{sample}" + FASTQ_EXT),
                    fastpHTML=join(PREPROC_DIR, "fastp.{sample}.html"),
                    fastpJSON=join(PREPROC_DIR, "fastp.{sample}.json"),
                log:
                    join(PREPROC_LOG_DIR, "fastqTrim.{sample}.log"),
                threads: THREADS
                run:
                    LOG_FQTRIM.info("Running fastqTrim.smk...")
                    startTime = time()

                    shell(
                    """
                    fastp \
                     --thread {threads} \
                     --in1 {input.fastq} \
                     --out1 {output.fastqClean} \
                     --html {output.fastpHTML} \
                     --json {output.fastpJSON} \
                     &> {log}
                    """
                    )

                    endTime = time()
                    elapseTime = endTime - startTime
                    LOG_FQTRIM.info(
                        f"Preprocessed input FASTQ file of {wildcards.sample} in {elapseTime:.2f} seconds"
                    )

        makedirs(PREPROC_DIR, exist_ok=True)
        df = pd.read_csv(BATCH, sep="\t", header=None)
        if df.shape[1] == 3:
            df.columns = ["cond", "rep", "fq"]
            df["fq"] = df["fq"].apply(
                lambda path: utils.process_path(path, PREPROC_DIR)
            )
        elif df.shape[1] == 4:
            df.columns = ["cond", "rep", "fq1", "fq2"]
            df["fq1"] = df["fq1"].apply(
                lambda path: utils.process_path(path, PREPROC_DIR)
            )
            df["fq2"] = df["fq2"].apply(
                lambda path: utils.process_path(path, PREPROC_DIR)
            )
        else:
            raise ValueError(f"Unexpected number of columns in {BATCH}")
        df.to_csv(BATCH_UDT, sep="\t", index=False, header=False)

    else:
        if PAIREDEND:

            rule fastp_cleanFQ:
                """
                Preprocess input FASTQ file
                """
                input:
                    fastqLeft=FASTQ_LEFT,
                    fastqRight=FASTQ_RIGHT,
                output:
                    fastqCleanLeft=join(PREPROC_DIR, SAMPLE + READS[0] + FASTQ_EXT),
                    fastqCleanRight=join(PREPROC_DIR, SAMPLE + READS[1] + FASTQ_EXT),
                    fastpHTML=join(PREPROC_DIR, "fastp.html"),
                    fastpJSON=join(PREPROC_DIR, "fastp.json"),
                log:
                    join(PREPROC_LOG_DIR, "fastqTrim.log"),
                threads: THREADS
                run:
                    LOG_FQTRIM.info("Running fastqTrim.smk...")
                    startTime = time()

                    shell(
                    """
                    fastp \
                     --thread {threads} \
                     --in1 {input.fastqLeft} --in2 {input.fastqRight} \
                     --out1 {output.fastqCleanLeft} --out2 {output.fastqCleanRight} \
                     --html {output.fastpHTML} \
                     --json {output.fastpJSON} \
                     &> {log}
                    """
                    )

                    endTime = time()
                    elapseTime = endTime - startTime
                    LOG_FQTRIM.info(
                        f"Preprocessed input FASTQ file in {elapseTime:.2f} seconds"
                    )

        else:

            rule fastp_cleanFQ:
                """
                Preprocess input FASTQ file
                """
                input:
                    fastq=FASTQ,
                output:
                    fastqClean=join(PREPROC_DIR, SAMPLE + FASTQ_EXT),
                    fastpHTML=join(PREPROC_DIR, "fastp.html"),
                    fastpJSON=join(PREPROC_DIR, "fastp.json"),
                log:
                    join(PREPROC_LOG_DIR, "fastqTrim.log"),
                threads: THREADS
                run:
                    LOG_FQTRIM.info("Running fastqTrim.smk...")
                    startTime = time()

                    shell(
                    """
                    fastp \
                     --thread {threads} \
                     --in1 {input.fastq} \
                     --out1 {output.fastqClean} \
                     --html {output.fastpHTML} \
                     --json {output.fastpJSON} \
                     &> {log}
                    """
                    )

                    endTime = time()
                    elapseTime = endTime - startTime
                    LOG_FQTRIM.info(
                        f"Preprocessed input FASTQ file in {elapseTime:.2f} seconds"
                    )

else:
    if BATCH:
        if PAIREDEND:

            rule copy_fastq:
                """
                Copy FASTQ file(s) to objective directory
                """
                input:
                    fastqLeft=lambda wildcards: SAMPLE_DICT[wildcards.sample][0],
                    fastqRight=lambda wildcards: SAMPLE_DICT[wildcards.sample][1],
                output:
                    fastqCleanLeft=join(PREPROC_DIR, "{sample}" + READS[0] + FASTQ_EXT),
                    fastqCleanRight=join(PREPROC_DIR, "{sample}" + READS[1] + FASTQ_EXT),
                log:
                    join(PREPROC_LOG_DIR, "fastqTrim.{sample}.log"),
                threads: THREADS
                run:
                    LOG_FQTRIM.info("Running fastqTrim.smk...")
                    startTime = time()

                    makedirs(PREPROC_DIR, exist_ok=True)

                    shell(
                    """
                    ln -sr {input.fastqLeft} {output.fastqCleanLeft} &> {log}
                    ln -sr {input.fastqRight} {output.fastqCleanRight} &>> {log}
                    """
                    )

                    endTime = time()
                    elapseTime = endTime - startTime
                    LOG_FQTRIM.info(
                        f"Copied FASTQ file(s) of {wildcards.sample} to objective directory {elapseTime:.2f} seconds"
                    )

        else:

            rule copy_fastq:
                """
                Copy FASTQ file(s) to objective directory
                """
                input:
                    fastq=lambda wildcards: SAMPLE_DICT[wildcards.sample][0],
                output:
                    fastqClean=join(PREPROC_DIR, "{sample}" + FASTQ_EXT),
                log:
                    join(PREPROC_LOG_DIR, "fastqTrim.{sample}.log"),
                threads: THREADS
                run:
                    LOG_FQTRIM.info("Running fastqTrim.smk...")
                    startTime = time()

                    makedirs(PREPROC_DIR, exist_ok=True)

                    shell(
                    """
                    ln -sr {input.fastq} {output.fastqClean} &> {log}
                    """
                    )

                    endTime = time()
                    elapseTime = endTime - startTime
                    LOG_FQTRIM.info(
                        f"Copied FASTQ file(s) of {wildcards.sample} to objective directory {elapseTime:.2f} seconds"
                    )

        makedirs(PREPROC_DIR, exist_ok=True)
        df = pd.read_csv(BATCH, sep="\t", header=None)
        if df.shape[1] == 3:
            df.columns = ["cond", "rep", "fq"]
            df["fq"] = df["fq"].apply(
                lambda path: utils.process_path(path, PREPROC_DIR)
            )
        elif df.shape[1] == 4:
            df.columns = ["cond", "rep", "fq1", "fq2"]
            df["fq1"] = df["fq1"].apply(
                lambda path: utils.process_path(path, PREPROC_DIR)
            )
            df["fq2"] = df["fq2"].apply(
                lambda path: utils.process_path(path, PREPROC_DIR)
            )
        else:
            raise ValueError(f"Unexpected number of columns in {BATCH}")
        df.to_csv(BATCH_UDT, sep="\t", index=False, header=False)

    else:
        if PAIREDEND:

            rule copy_fastq:
                """
                Copy FASTQ file(s) to objective directory
                """
                input:
                    fastqLeft=FASTQ_LEFT,
                    fastqRight=FASTQ_RIGHT,
                output:
                    fastqCleanLeft=join(PREPROC_DIR, SAMPLE + READS[0] + FASTQ_EXT),
                    fastqCleanRight=join(PREPROC_DIR, SAMPLE + READS[1] + FASTQ_EXT),
                log:
                    join(PREPROC_LOG_DIR, "fastqTrim.log"),
                threads: THREADS
                run:
                    LOG_FQTRIM.info("Running fastqTrim.smk...")
                    startTime = time()

                    makedirs(PREPROC_DIR, exist_ok=True)

                    shell(
                    """
                    ln -sr {input.fastqLeft} {output.fastqCleanLeft} &> {log}
                    ln -sr {input.fastqRight} {output.fastqCleanRight} &>> {log}
                    """
                    )

                    endTime = time()
                    elapseTime = endTime - startTime
                    LOG_FQTRIM.info(
                        f"Copied FASTQ file(s) to objective directory {elapseTime:.2f} seconds"
                    )

        else:

            rule copy_fastq:
                """
                Copy FASTQ file(s) to objective directory
                """
                input:
                    fastq=FASTQ,
                output:
                    fastqClean=join(PREPROC_DIR, SAMPLE + FASTQ_EXT),
                log:
                    join(PREPROC_LOG_DIR, "fastqTrim.log"),
                threads: THREADS
                run:
                    LOG_FQTRIM.info("Running fastqTrim.smk...")

                    makedirs(PREPROC_DIR, exist_ok=True)

                    shell(
                    """
                    ln -sr {input.fastq} {output.fastqClean} &> {log}
                    """
                    )

                    endTime = time()
                    elapseTime = endTime - startTime
                    LOG_FQTRIM.info(
                        f"Copied FASTQ file(s) to objective directory {elapseTime:.2f} seconds"
                    )
