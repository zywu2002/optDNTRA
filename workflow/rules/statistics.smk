from os.path import join
import logging
import common_functions as cf
from time import time
from prettytable import PrettyTable


### Logging Setup ### ----------------------------------------------
LOG_STATS = cf.get_logger("STATISTICS", VERBOSE)


### Rule ### -------------------------------------------------------
rule statistics:
    """
    Generate statistics for Optimization of De Novo Transcriptome Assembly (OptDNTA)
    """
    input:
        transcript = join(PREPROC_DIR, "transcript.fa"),
        transcriptFlt1 = join(RMREDUNDANT_DIR, "transcript.flt1.fa"),
        transcriptFlt2 = join(ORFPREDICT_DIR, "transcript.flt2.fa"),
        transcriptFlt3 = join(TRANSEVID_DIR, "transcript.flt.final.fa")
    output:
        log = join(STATS_LOG_DIR, "statistics.log")
    run:
        LOG_STATS.info("Running statistics.smk...")
        startTime = time()

        seqNumFlt0 = str(cf.count_fastx(input.transcript, LOG_STATS, identifier=">"))
        seqNumFlt1 = str(cf.count_fastx(input.transcriptFlt1, LOG_STATS, identifier=">"))
        seqNumFlt2 = str(cf.count_fastx(input.transcriptFlt2, LOG_STATS, identifier=">"))
        seqNumFlt3 = str(cf.count_fastx(input.transcriptFlt3, LOG_STATS, identifier=">"))

        with open(output.log, "w") as FOUT:
            FOUT.write("Transcript" + "\t" + "Num" + "\n")
            FOUT.write("transcript.fa" + "\t" + seqNumFlt0 + "\n")
            FOUT.write("transcript.flt1.fa" + "\t" + seqNumFlt1 + "\n")
            FOUT.write("transcript.flt2.fa" + "\t" + seqNumFlt2 + "\n")
            FOUT.write("transcript.flt.final.fa" + "\t" + seqNumFlt3)

        table = PrettyTable()
        table.field_names = ["Transcript", "Num"]
        table.add_row(["transcript.fa", seqNumFlt0])
        table.add_row(["transcript.flt1.fa", seqNumFlt1])
        table.add_row(["transcript.flt2.fa", seqNumFlt2])
        table.add_row(["transcript.flt.final.fa", seqNumFlt3])

        table.align = "l"
        table.header = True
        table.border = True

        LOG_STATS.info("\n" + str(table))

        endTime = time()
        elapseTime = endTime - startTime
        LOG_STATS.info(f"Generated statistics in {elapseTime:.2f} seconds")