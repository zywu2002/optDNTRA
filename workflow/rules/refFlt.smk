from os.path import join
import logging
import common_functions as cf
from time import time


### Configuration Setup ### ----------------------------------------
REFCDNA = config["reference"]


### Logging Setup ### ----------------------------------------
LOG_REFFLT = cf.get_logger("REFFLT", VERBOSE)


### Rule ### ----------------------------------------
if REFERENCE:
    rule filter_based_on_reference:
        """
        Perform BLAST search against reference cDNA and filter transcripts
        """
        input:
            transcript = join(ORFPREDICT_DIR, "transcript.flt2.fa"),
            refcDNA = REFCDNA
        output:
            transcriptFlt = temp(join(TRANSEVID_DIR, "transcript.noBlast1.fasta")),
            refBlastOut = join(TRANSEVID_DIR, "refBlastn.outfmt6"),
            refBlastID = temp(join(TRANSEVID_DIR, "blast.ref_cdna.id.txt")),
            transcriptFltID = temp(join(TRANSEVID_DIR, "transcriptFlt.id.txt")),
            transcriptFltNoBlastID = temp(join(TRANSEVID_DIR, "transcriptFlt.noBlast.id.txt")),
        log:
            join(OPTIM_LOG_DIR, "refFlt.log")
        threads:
            THREADS
        run:
            LOG_REFFLT.info("Running refFlt.smk...")
            startTime = time()

            shell("""
            makeblastdb -in {input.refcDNA} -dbtype nucl &> {log}
            
            blastn \
            -query {input.transcript} \
            -db {input.refcDNA} \
            -out {output.refBlastOut} \
            -evalue 1e-5 \
            -max_target_seqs 1 \
            -outfmt 6 \
            -num_threads {threads} \
            &>> {log}

            cut -f1 {output.refBlastOut} | sort | uniq > {output.refBlastID} 2>> {log}
            awk '/^>/ {{sub("^>", ""); sub(" .*", ""); print}}' {input.transcript} | sort > {output.transcriptFltID} 2>> {log}
            comm -13 {output.refBlastID} {output.transcriptFltID} | sort > {output.transcriptFltNoBlastID} 2>> {log}
            seqtk subseq {input.transcript} {output.transcriptFltNoBlastID} > {output.transcriptFlt} 2>> {log}        
            """)

            endTime = time()
            elapseTime = endTime - startTime
            LOG_REFFLT.info(f"Performed BLAST search against reference cDNA and filtered transcripts in {elapseTime:.2f} seconds")
else:
    rule filter_based_on_reference:
        """
        Perform BLAST search against reference cDNA and filter transcripts
        """
        input:
            transcript = join(ORFPREDICT_DIR, "transcript.flt2.fa"),
        output:
            transcriptFlt = temp(join(TRANSEVID_DIR, "transcript.noBlast1.fasta")),
            transcriptFltID = temp(join(TRANSEVID_DIR, "transcriptFlt.id.txt")),
            transcriptFltNoBlastID = temp(join(TRANSEVID_DIR, "transcriptFlt.noBlast.id.txt"))
        log:
            join(OPTIM_LOG_DIR, "refFlt.log")
        threads:
            THREADS
        run:
            shell("""
            awk '/^>/ {{sub("^>", ""); sub(" .*", ""); print}}' {input.transcript} | sort > {output.transcriptFltID} 2> {log}
            cp {output.transcriptFltID} {output.transcriptFltNoBlastID} 2>> {log}
            cp {input.transcript} {output.transcriptFlt} 2>> {log}
            """)