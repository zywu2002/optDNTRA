from os.path import join
from os import makedirs
import logging
import common_functions as cf
from time import time
import pandas as pd
import shutil


### Logging Setup ### ----------------------------------------
LOG_EXPRFLT = cf.get_logger("EXPRFLT", VERBOSE)


### Rule ### ----------------------------------------
if BATCH:

    checkpoint estimate_tx_abundance:
        """
        Estimating transcript abundance
        """
        input:
            batch=BATCH_UDT if TRIM else BATCH,
            transcript=join(TRANSEVID_DIR, "transcript.noBlast2.fasta"),
        output:
            complete=temp(join(TRANSEVID_DIR, "estimate_tx_abundance.complete")),
        log:
            join(OPTIM_LOG_DIR, "exprFlt-estimate_tx_abundance.log"),
        params:
            trinityToolkitOutdir=join(TRANSEVID_DIR, "trinityToolkitOut"),
            ssLibType=SSLIBTYPE,
        threads: THREADS
        run:
            LOG_EXPRFLT.info("Running exprFlt.smk...")
            startTime = time()

            makedirs(params.trinityToolkitOutdir, exist_ok=True)

            shell(
            """
            align_and_estimate_abundance.pl \
            --transcripts {input.transcript} \
            --seqType fq \
            --samples_file {input.batch} \
            --est_method RSEM \
            --aln_method bowtie2 \
            --SS_lib_type {params.ssLibType} \
            --prep_reference \
            --thread_count {threads} \
            &> {log}
            """
            )

            df = pd.read_csv(input.batch, sep="\t", header=None)
            rsemResDirs = [dir for dir in df.iloc[:, 1].tolist()]
            for dir in rsemResDirs:
                shutil.move(dir, join(params.trinityToolkitOutdir, dir))

            shell(
            """
            touch {output.complete}
            """
            )

            endTime = time()
            elapseTime = endTime - startTime
            LOG_ORFPRED.info(
                f"Estimated transcript abundance in {elapseTime:.2f} seconds."
            )

    def get_RSEMIsoRes_batch(wildcards):
        checkpointOut = checkpoints.estimate_tx_abundance.get(
            **wildcards
        ).output.complete
        batchFile = BATCH_UDT if TRIM else BATCH
        df = pd.read_csv(batchFile, sep="\t", header=None)
        return [
            join(TRANSEVID_DIR, "trinityToolkitOut", dir, "RSEM.isoforms.results")
            for dir in df.iloc[:, 1].tolist()
        ]

    checkpoint build_expr_matrics:
        """
        Build transcript and gene expression matrics
        """
        input:
            RSEMIsoRes=get_RSEMIsoRes_batch,
        output:
            RSEMIsoTPM=join(
                TRANSEVID_DIR, "trinityToolkitOut", "RSEM.isoform.TMM.EXPR.matrix"
            ),
        log:
            join(OPTIM_LOG_DIR, "exprFlt-build_expr_matrics.log"),
        params:
            RSEMprefix=join(TRANSEVID_DIR, "trinityToolkitOut", "RSEM"),
        threads: THREADS
        run:
            startTime = time()

            RSEMIsoRes = " ".join(input.RSEMIsoRes)

            shell(
            """
            abundance_estimates_to_matrix.pl \
            --est_method RSEM \
            --gene_trans_map none \
            --name_sample_by_basedir \
            --out_prefix {params.RSEMprefix} \
            {RSEMIsoRes} \
            &> {log}

            rm {TRANSEVID_DIR}/estimate_tx_abundance.complete 2>> {log}
            """
            )

            endTime = time()
            elapseTime = endTime - startTime
            LOG_ORFPRED.info(
                f"Built transcript and gene expression matrics in {elapseTime:.2f} seconds."
            )

    def get_RSEMIsoTPM_batch(wildcards):
        checkpointOut = checkpoints.build_expr_matrics.get(
            **wildcards
        ).output.RSEMIsoTPM
        return checkpointOut

    rule filter_based_on_exprLevel:
        """
        Filtering transcripts based on expression values
        """
        input:
            transcript1=join(TRANSEVID_DIR, "transcript.noBlast2.fasta"),
            transcript2=join(ORFPREDICT_DIR, "transcript.flt2.fa"),
            RSEMIsoTPM=get_RSEMIsoTPM_batch,
            transcriptFltNoBlast2ID=join(TRANSEVID_DIR, "transcriptFlt.noBlast2.id.txt"),
            transcriptFltID=join(TRANSEVID_DIR, "transcriptFlt.id.txt"),
            transcriptPep=join(
                ORFPREDICT_DIR, "transDecoder", "transcript.flt1.fa.transdecoder.pep"
            ),
        output:
            transcriptNoBlast2TPM1fa=temp(
                join(TRANSEVID_DIR, "transcript.noBlast2.tpm1.fasta")
            ),
            transcriptNoBlast2TPM1ID=temp(
                join(TRANSEVID_DIR, "transcript.noBlast2.tpm1.id.txt")
            ),
            transcriptNoBlast2NoTPM1ID=temp(
                join(TRANSEVID_DIR, "transcript.noBlast2.notpm1.id.txt")
            ),
            transcriptFltFinalID=temp(join(TRANSEVID_DIR, "transcriptFltFinal.id.txt")),
            transcriptFlt1Pep=temp(
                join(ORFPREDICT_DIR, "transDecoder", "transcript.flt1.pep")
            ),
            transcriptFlt=join(TRANSEVID_DIR, "transcript.flt.final.fa"),
            transcriptFltPep=join(TRANSEVID_DIR, "transcript.flt.final.pep"),
        log:
            join(OPTIM_LOG_DIR, "exprFlt-filter_based_on_exprLevel.log"),
        threads: THREADS
        run:
            startTime = time()

            shell(
            """
            filter_low_expr_transcripts.pl \
            --matrix {input.RSEMIsoTPM} \
            --transcripts {input.transcript1} \
            --min_expr_any 1 \
            > {output.transcriptNoBlast2TPM1fa} \
            2> {log}
            """
            )

            shell(
                """
                        awk '/^>/ {{sub("^>", ""); sub(" .*", ""); print}}' {output.transcriptNoBlast2TPM1fa} | sort > {output.transcriptNoBlast2TPM1ID} 2>> {log}
                        comm -13 {output.transcriptNoBlast2TPM1ID} {input.transcriptFltNoBlast2ID} | sort > {output.transcriptNoBlast2NoTPM1ID} 2>> {log}
                        comm -13 {output.transcriptNoBlast2NoTPM1ID} {input.transcriptFltID} | sort > {output.transcriptFltFinalID} 2>> {log}
                        seqtk subseq {input.transcript2} {output.transcriptFltFinalID} > {output.transcriptFlt} 2>> {log}

                        sed -E 's/^>([^.]+)\\.p.*/>\\1/' {input.transcriptPep} > {output.transcriptFlt1Pep} 2>> {log}
                        seqtk subseq {output.transcriptFlt1Pep} {output.transcriptFltFinalID} > {output.transcriptFltPep} 2>> {log}        

                        rm {TRANSEVID_DIR}/transcript.noBlast2.fasta.* 2>> {log} 
                        """
            )

            endTime = time()
            elapseTime = endTime - startTime
            LOG_EXPRFLT.info(
                f"Filtered transcripts based on expression values in {elapseTime:.2f} seconds."
            )

else:
    if PAIREDEND:

        checkpoint estimate_tx_abundance:
            """
            Estimating transcript abundance
            """
            input:
                fqLeft=FASTQ_LEFT,
                fqRight=FASTQ_RIGHT,
                transcript=join(TRANSEVID_DIR, "transcript.noBlast2.fasta"),
            output:
                RSEMIsoRes=join(
                    TRANSEVID_DIR, "trinityToolkitOut", "RSEM.isoforms.results"
                ),
            log:
                join(OPTIM_LOG_DIR, "exprFlt-estimate_tx_abundance.log"),
            params:
                trinityToolkitOutdir=join(TRANSEVID_DIR, "trinityToolkitOut"),
                ssLibType=SSLIBTYPE,
            threads: THREADS
            run:
                LOG_EXPRFLT.info("Running exprFlt.smk...")
                startTime = time()

                makedirs(params.trinityToolkitOutdir, exist_ok=True)

                shell(
                """
                align_and_estimate_abundance.pl \
                --transcripts {input.transcript} \
                --seqType fq \
                --left {input.fqLeft} \
                --right {input.fqRight} \
                --est_method RSEM \
                --aln_method bowtie2 \
                --SS_lib_type {params.ssLibType} \
                --prep_reference \
                --output_dir {params.trinityToolkitOutdir} \
                --thread_count {threads} \
                &> {log}
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_ORFPRED.info(
                    f"Estimated transcript abundance in {elapseTime:.2f} seconds."
                )

        def get_RSEMIsoRes(wildcards):
            checkpointOut = checkpoints.estimate_tx_abundance.get(
                **wildcards
            ).output.RSEMIsoRes
            return checkpointOut

        checkpoint build_expr_matrics:
            """
            Build transcript and gene expression matrics
            """
            input:
                RSEMIsoRes=get_RSEMIsoRes,
            output:
                RSEMIsoTPM=join(
                    TRANSEVID_DIR,
                    "trinityToolkitOut",
                    "RSEM.isoform.TPM.not_cross_norm",
                ),
            log:
                join(OPTIM_LOG_DIR, "exprFlt-build_expr_matrics.log"),
            params:
                RSEMprefix=join(TRANSEVID_DIR, "trinityToolkitOut", "RSEM"),
            threads: THREADS
            run:
                startTime = time()

                shell(
                """
                abundance_estimates_to_matrix.pl \
                --est_method RSEM \
                --gene_trans_map none \
                --out_prefix {params.RSEMprefix} \
                {input.RSEMIsoRes} \
                &> {log}
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_ORFPRED.info(
                    f"Built transcript and gene expression matrics in {elapseTime:.2f} seconds."
                )

        def get_RSEMIsoTPM(wildcards):
            checkpointOut = checkpoints.build_expr_matrics.get(
                **wildcards
            ).output.RSEMIsoTPM
            return checkpointOut

        rule filter_based_on_exprLevel:
            """
            Filtering transcripts based on expression values
            """
            input:
                transcript1=join(TRANSEVID_DIR, "transcript.noBlast2.fasta"),
                transcript2=join(ORFPREDICT_DIR, "transcript.flt2.fa"),
                RSEMIsoTPM=get_RSEMIsoTPM,
                transcriptFltNoBlast2ID=join(
                    TRANSEVID_DIR, "transcriptFlt.noBlast2.id.txt"
                ),
                transcriptFltID=join(TRANSEVID_DIR, "transcriptFlt.id.txt"),
                transcriptPep=join(
                    ORFPREDICT_DIR,
                    "transDecoder",
                    "transcript.flt1.fa.transdecoder.pep",
                ),
            output:
                transcriptNoBlast2TPM1fa=temp(
                    join(TRANSEVID_DIR, "transcript.noBlast2.tpm1.fasta")
                ),
                transcriptNoBlast2TPM1ID=temp(
                    join(TRANSEVID_DIR, "transcript.noBlast2.tpm1.id.txt")
                ),
                transcriptNoBlast2NoTPM1ID=temp(
                    join(TRANSEVID_DIR, "transcript.noBlast2.notpm1.id.txt")
                ),
                transcriptFltFinalID=temp(
                    join(TRANSEVID_DIR, "transcriptFltFinal.id.txt")
                ),
                transcriptFlt1Pep=temp(
                    join(ORFPREDICT_DIR, "transDecoder", "transcript.flt1.pep")
                ),
                transcriptFlt=join(TRANSEVID_DIR, "transcript.flt.final.fa"),
                transcriptFltPep=join(TRANSEVID_DIR, "transcript.flt.final.pep"),
            log:
                join(OPTIM_LOG_DIR, "exprFlt-filter_based_on_exprLevel.log"),
            threads: THREADS
            run:
                startTime = time()

                shell(
                """
                filter_low_expr_transcripts.pl \
                --matrix {input.RSEMIsoTPM} \
                --transcripts {input.transcript1} \
                --min_expr_any 1 \
                > {output.transcriptNoBlast2TPM1fa} \
                2> {log}
                """
                )

                shell(
                """
                awk '/^>/ {{sub("^>", ""); sub(" .*", ""); print}}' {output.transcriptNoBlast2TPM1fa} | sort > {output.transcriptNoBlast2TPM1ID} 2>> {log}
                comm -13 {output.transcriptNoBlast2TPM1ID} {input.transcriptFltNoBlast2ID} | sort > {output.transcriptNoBlast2NoTPM1ID} 2>> {log}
                comm -13 {output.transcriptNoBlast2NoTPM1ID} {input.transcriptFltID} | sort > {output.transcriptFltFinalID} 2>> {log}
                seqtk subseq {input.transcript2} {output.transcriptFltFinalID} > {output.transcriptFlt} 2>> {log}

                sed -E 's/^>([^.]+)\\.p.*/>\\1/' {input.transcriptPep} > {output.transcriptFlt1Pep} 2>> {log}
                seqtk subseq {output.transcriptFlt1Pep} {output.transcriptFltFinalID} > {output.transcriptFltPep} 2>> {log}        

                rm {TRANSEVID_DIR}/transcript.noBlast2.fasta.* 2>> {log} 
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_EXPRFLT.info(
                    f"Filtered transcripts based on expression values in {elapseTime:.2f} seconds."
                )

    else:

        checkpoint estimate_tx_abundance:
            """
            Estimating transcript abundance
            """
            input:
                fastq=FASTQ,
                transcript=join(TRANSEVID_DIR, "transcript.noBlast2.fasta"),
            output:
                RSEMIsoRes=join(
                    TRANSEVID_DIR, "trinityToolkitOut", "RSEM.isoforms.results"
                ),
            log:
                join(OPTIM_LOG_DIR, "exprFlt-estimate_tx_abundance.log"),
            params:
                trinityToolkitOutdir=join(TRANSEVID_DIR, "trinityToolkitOut"),
                ssLibType=SSLIBTYPE,
            threads: THREADS
            run:
                LOG_EXPRFLT.info("Running exprFlt.smk...")
                startTime = time()

                makedirs(params.trinityToolkitOutdir, exist_ok=True)

                shell(
                """
                align_and_estimate_abundance.pl \
                --transcripts {input.transcript} \
                --seqType fq \
                --single {input.fastq} \
                --est_method RSEM \
                --aln_method bowtie2 \
                --SS_lib_type {params.ssLibType} \
                --prep_reference \
                --output_dir {params.trinityToolkitOutdir} \
                --thread_count {threads} \
                &> {log}
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_ORFPRED.info(
                    f"Estimated transcript abundance in {elapseTime:.2f} seconds."
                )

        def get_RSEMIsoRes(wildcards):
            checkpointOut = checkpoints.estimate_tx_abundance.get(
                **wildcards
            ).output.RSEMIsoRes
            return checkpointOut

        checkpoint build_expr_matrics:
            """
            Build transcript and gene expression matrics
            """
            input:
                RSEMIsoRes=get_RSEMIsoRes,
            output:
                RSEMIsoTPM=join(
                    TRANSEVID_DIR,
                    "trinityToolkitOut",
                    "RSEM.isoform.TPM.not_cross_norm",
                ),
            log:
                join(OPTIM_LOG_DIR, "exprFlt-build_expr_matrics.log"),
            params:
                RSEMprefix=join(TRANSEVID_DIR, "trinityToolkitOut", "RSEM"),
            threads: THREADS
            run:
                startTime = time()

                shell(
                """
                abundance_estimates_to_matrix.pl \
                --est_method RSEM \
                --gene_trans_map none \
                --out_prefix {params.RSEMprefix} \
                {input.RSEMIsoRes} \
                &> {log}
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_ORFPRED.info(
                    f"Built transcript and gene expression matrics in {elapseTime:.2f} seconds."
                )

        def get_RSEMIsoTPM(wildcards):
            checkpointOut = checkpoints.build_expr_matrics.get(
                **wildcards
            ).output.RSEMIsoTPM
            return checkpointOut

        rule filter_based_on_exprLevel:
            """
            Filtering transcripts based on expression values
            """
            input:
                transcript1=join(TRANSEVID_DIR, "transcript.noBlast2.fasta"),
                transcript2=join(ORFPREDICT_DIR, "transcript.flt2.fa"),
                RSEMIsoTPM=get_RSEMIsoTPM,
                transcriptFltNoBlast2ID=join(
                    TRANSEVID_DIR, "transcriptFlt.noBlast2.id.txt"
                ),
                transcriptFltID=join(TRANSEVID_DIR, "transcriptFlt.id.txt"),
                transcriptPep=join(
                    ORFPREDICT_DIR,
                    "transDecoder",
                    "transcript.flt1.fa.transdecoder.pep",
                ),
            output:
                transcriptNoBlast2TPM1fa=temp(
                    join(TRANSEVID_DIR, "transcript.noBlast2.tpm1.fasta")
                ),
                transcriptNoBlast2TPM1ID=temp(
                    join(TRANSEVID_DIR, "transcript.noBlast2.tpm1.id.txt")
                ),
                transcriptNoBlast2NoTPM1ID=temp(
                    join(TRANSEVID_DIR, "transcript.noBlast2.notpm1.id.txt")
                ),
                transcriptFltFinalID=temp(
                    join(TRANSEVID_DIR, "transcriptFltFinal.id.txt")
                ),
                transcriptFlt1Pep=temp(
                    join(ORFPREDICT_DIR, "transDecoder", "transcript.flt1.pep")
                ),
                transcriptFlt=join(TRANSEVID_DIR, "transcript.flt.final.fa"),
                transcriptFltPep=join(TRANSEVID_DIR, "transcript.flt.final.pep"),
            log:
                join(OPTIM_LOG_DIR, "exprFlt-filter_based_on_exprLevel.log"),
            threads: THREADS
            run:
                startTime = time()

                shell(
                """
                filter_low_expr_transcripts.pl \
                --matrix {input.RSEMIsoTPM} \
                --transcripts {input.transcript1} \
                --min_expr_any 1 \
                > {output.transcriptNoBlast2TPM1fa} \
                2> {log}
                """
                )

                shell(
                """
                awk '/^>/ {{sub("^>", ""); sub(" .*", ""); print}}' {output.transcriptNoBlast2TPM1fa} | sort > {output.transcriptNoBlast2TPM1ID} 2>> {log}
                comm -13 {output.transcriptNoBlast2TPM1ID} {input.transcriptFltNoBlast2ID} | sort > {output.transcriptNoBlast2NoTPM1ID} 2>> {log}
                comm -13 {output.transcriptNoBlast2NoTPM1ID} {input.transcriptFltID} | sort > {output.transcriptFltFinalID} 2>> {log}
                seqtk subseq {input.transcript2} {output.transcriptFltFinalID} > {output.transcriptFlt} 2>> {log}

                sed -E 's/^>([^.]+)\\.p.*/>\\1/' {input.transcriptPep} > {output.transcriptFlt1Pep} 2>> {log}
                seqtk subseq {output.transcriptFlt1Pep} {output.transcriptFltFinalID} > {output.transcriptFltPep} 2>> {log}        

                rm {TRANSEVID_DIR}/transcript.noBlast2.fasta.* 2>> {log} 
                """
                )

                endTime = time()
                elapseTime = endTime - startTime
                LOG_EXPRFLT.info(
                    f"Filtered transcripts based on expression values in {elapseTime:.2f} seconds."
                )
