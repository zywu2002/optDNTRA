from os.path import join
import logging
import common_functions as cf
from time import time


### Configuration Setup ### ----------------------------------------

# BlastP Search
SWISSPROT = config["swiss_prot"]

# Pfam search
PFAMHMM = config["pfam_hmm"]


### Logging Setup ### ----------------------------------------
LOG_ORFPRED = cf.get_logger("ORFPRED", VERBOSE)


### Rule ### ----------------------------------------
if SSLIBTYPE:
    checkpoint transdecoder_longorfs:
        """
        Extract long open reading frames (ORFs) using TransDecoder.LongOrfs
        """
        input:
            transcript = join(RMREDUNDANT_DIR, "transcript.flt1.fa")
        output:
            longestPep = join(ORFPREDICT_DIR, "transDecoder", "transcript.flt1.fa.transdecoder_dir", "longest_orfs.pep")
        log:
            join(OPTIM_LOG_DIR, "orfPred-transdecoder_longorfs.log")
        params:
            outDir = join(ORFPREDICT_DIR, "transDecoder"),
        threads:
            THREADS
        run:
            LOG_ORFPRED.info("Running orfPred.smk...")
            startTime = time()

            shell("""
            TransDecoder.LongOrfs \
            -t {input.transcript} \
            -S \
            --output_dir {params.outDir} \
            &> {log}
            """)

            endTime = time()
            elapseTime = endTime - startTime
            LOG_ORFPRED.info(f"Extracted long open reading frames (ORFs) using TransDecoder.LongOrfs in {elapseTime:.2f} seconds")
else:
    checkpoint transdecoder_longorfs:
        """
        Extract long open reading frames (ORFs) using TransDecoder.LongOrfs
        """
        input:
            transcript = join(RMREDUNDANT_DIR, "transcript.flt1.fa")
        output:
            longestPep = join(ORFPREDICT_DIR, "transDecoder", "transcript.flt1.fa.transdecoder_dir", "longest_orfs.pep")
        log:
            join(OPTIM_LOG_DIR, "orfPred-transdecoder_longorfs.log")
        params:
            outDir = join(ORFPREDICT_DIR, "transDecoder"),
        threads:
            THREADS
        run:
            LOG_ORFPRED.info("Running orfPred.smk...")
            startTime = time()

            shell("""
            TransDecoder.LongOrfs \
            -t {input.transcript} \
            --output_dir {params.outDir} \
            &> {log}
            """)

            endTime = time()
            elapseTime = endTime - startTime
            LOG_ORFPRED.info(f"Extracted long open reading frames (ORFs) using TransDecoder.LongOrfs in {elapseTime:.2f} seconds")

def get_longestORFs_pep(wildcards):
    checkpointOut = checkpoints.transdecoder_longorfs.get(**wildcards).output.longestPep
    return checkpointOut

rule blastp_search:
    """
    Perform BlastP search using Diamond blastp
    """
    input:
        swissProt = SWISSPROT,
        longestOrfs = get_longestORFs_pep
    output:
        blastpOut = join(ORFPREDICT_DIR, "transcript.flt1.blastp.outfmt6")
    log:
        join(OPTIM_LOG_DIR, "orfPred-blastp_search.log")
    threads:
        THREADS
    run:
        startTime = time()

        shell("""
        diamond makedb \
         --in {input.swissProt} \
         --db {input.swissProt} \
         --threads {threads} \
         &> {log}

        diamond blastp \
         --db {input.swissProt} \
         --query {input.longestOrfs} \
         --out {output.blastpOut} \
         --evalue 1e-5 \
         --max-target-seqs 1 \
         --outfmt 6 \
         --threads {threads} \
         &>> {log}
        """)

        endTime = time()
        elapseTime = endTime - startTime
        LOG_ORFPRED.info(f"Performed BlastP search using Diamond blastp in {elapseTime:.2f} seconds")

rule pfam_search:
    """
    Perform Pfam domain search using hmmsearch
    """
    input:
        pfamHMM = PFAMHMM,
        longestOrfs = get_longestORFs_pep
    output:
        pfamOut = join(ORFPREDICT_DIR, "transcript.flt1.pfam.domtblout")
    log:
        join(OPTIM_LOG_DIR, "orfPred-pfam_search.log")
    threads:
        THREADS
    run:
        startTime = time()

        shell("""
        hmmsearch \
         --cpu {threads} \
         --domtblout {output.pfamOut} \
         -E 1e-10 \
         {input.pfamHMM} \
         {input.longestOrfs} \
         &> {log}
        """)

        endTime = time()
        elapseTime = endTime - startTime
        LOG_ORFPRED.info(f"Performed Pfam domain search using hmmsearch in {elapseTime:.2f} seconds")

checkpoint transdecoder_predict:
    """
    Predict likely coding regions using TransDecoder.Predict.
    """
    input:
        transcript = join(RMREDUNDANT_DIR, "transcript.flt1.fa"),
        blastpOut = rules.blastp_search.output.blastpOut,
        pfamOut = rules.pfam_search.output.pfamOut
    output:
        predictPep = join(ORFPREDICT_DIR, "transDecoder", "transcript.flt1.fa.transdecoder.pep")
    log:
        join(OPTIM_LOG_DIR, "orfPred-transdecoder_predict.log")
    params:
        outDir = join(ORFPREDICT_DIR, "transDecoder")
    threads:
        THREADS
    run:
        startTime = time()

        shell("""
        TransDecoder.Predict \
         -t {input.transcript} \
         --output_dir {params.outDir} \
         --retain_pfam_hits {input.pfamOut} \
         --retain_blastp_hits {input.blastpOut} \
         --single_best_only \
         &> {log}
        """)

        endTime = time()
        elapseTime = endTime - startTime
        LOG_ORFPRED.info(f"Predicted likely coding regions using TransDecoder.Predict in {elapseTime:.2f} seconds")

def get_predict_pep(wildcards):
    checkpointOut = checkpoints.transdecoder_predict.get(**wildcards).output.predictPep
    return checkpointOut

rule filter_step:
    """
    Filter transcripts based on TransDecoder.Predict output
    """
    input:
        predictPEP = get_predict_pep,
        transcript = join(RMREDUNDANT_DIR, "transcript.flt1.fa")
    output:
        transcriptFlt = join(ORFPREDICT_DIR, "transcript.flt2.fa"),
        transdecoderID = temp(join(ORFPREDICT_DIR, "transdecoder_id.txt"))
    log:
        join(OPTIM_LOG_DIR, "orfPred-filter_step.log")
    threads:
        THREADS
    run:
        startTime = time()

        shell("""
        awk '/^>/ {{gsub(/\\..*/, "", $1); print $1}} !/^>/ {{print}}' {input.predictPEP} | \
        awk '/^>/ {{print substr($0, 2)}}' > {output.transdecoderID} 2> {log}
        
        seqtk subseq {input.transcript} {output.transdecoderID} > {output.transcriptFlt} 2>> {log}
        """)
        
        endTime = time()
        elapseTime = endTime - startTime
        LOG_ORFPRED.info(f"Filtered transcripts based on TransDecoder.Predict output in {elapseTime:.2f} seconds")