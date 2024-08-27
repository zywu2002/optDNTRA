from os.path import join
import logging
import common_functions as cf
from time import time


### Configuration Setup ### ----------------------------------------
SWISSPROT = config["swiss_prot"]


### Logging Setup ### ----------------------------------------
LOG_DBFLT = cf.get_logger("DBFLT", VERBOSE)


### Rule ### ----------------------------------------
rule filter_based_on_database:
    """
    Perform BLAST search against Swiss-Prot database and filter transcripts
    """
    input:
        swissProt = SWISSPROT,
        transcript = join(TRANSEVID_DIR, "transcript.noBlast1.fasta"),
        transcriptFltNoBlastID = join(TRANSEVID_DIR, "transcriptFlt.noBlast.id.txt")
    output:
        swissprotBlastxOut = join(TRANSEVID_DIR, "swissprotBlastx.outfmt6"),
        transcriptFlt = temp(join(TRANSEVID_DIR, "transcript.noBlast2.fasta")),
        swissprotBlastID = temp(join(TRANSEVID_DIR, "blast.swissprot.id.txt")),
        transcriptFltNoBlast2ID = temp(join(TRANSEVID_DIR, "transcriptFlt.noBlast2.id.txt"))
    log:
        join(OPTIM_LOG_DIR, "dbFlt.log")
    threads:
        THREADS
    run:
        LOG_DBFLT.info("Running dbFlt.smk...")
        startTime = time()

        shell("""
        diamond blastx \
         --db {input.swissProt} \
         --query {input.transcript} \
         --out {output.swissprotBlastxOut} \
         --evalue 1e-5 \
         --max-target-seqs 1 \
         --outfmt 6 \
         --threads {threads} \
         &> {log}

        cut -f1 {output.swissprotBlastxOut} | sort | uniq > {output.swissprotBlastID} 2>> {log}
        comm -13 {output.swissprotBlastID} {input.transcriptFltNoBlastID} | sort > {output.transcriptFltNoBlast2ID} 2>> {log}
        seqtk subseq {input.transcript} {output.transcriptFltNoBlast2ID} > {output.transcriptFlt} 2>> {log}
        """)

        endTime = time()
        elapseTime = endTime - startTime
        LOG_DBFLT.info(f"Performed BLAST search against Swiss-Prot database and filtered transcripts in {elapseTime:.2f} seconds")