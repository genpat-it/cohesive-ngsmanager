nextflow.enable.dsl=2

include { step_0SQ_rawreads__fastq } from '../steps/step_0SQ_rawreads__fastq'
include { step_1PP_trimming__trimmomatic } from '../steps/step_1PP_trimming__trimmomatic'
include { step_1PP_trimming__fastp } from '../steps/step_1PP_trimming__fastp'
include { step_3TX_class__kraken } from '../steps/step_3TX_class__kraken'
include { extractKey } from '../functions/common.nf'
include { getInput;hasFastqData;hasEnoughFastqData;isIlluminaPaired;isIonTorrent;isNanopore; } from '../functions/parameters.nf'
include { isBacterium } from '../functions/samplesheet.nf'

workflow module_reads_processing {
    take: 
      rawReads
    main:
        rawReads.branch {
            with_data: hasFastqData(it[1])
            no_reads: true
        }
        .set { rawreads_branched }
        step_0SQ_rawreads__fastq(rawreads_branched.with_data)        

        rawreads_branched.with_data.branch {
            illumina: isIlluminaPaired(it[1])
            ion: isIonTorrent(it[1])
            nanopore: isNanopore(it[1])
            other: true // won't be processed
        }
        .set { trimming_by_seqtype }

        trimming_by_seqtype.illumina.branch {
            bacteria: isBacterium(it)
            other: true 
        }
        .set { trimming_illumina }

        // trimmomatic
        trimmed_by_trimmomatic = step_1PP_trimming__trimmomatic(trimming_illumina.other).trimmed

        // fastp
        trimmed_by_fastp = step_1PP_trimming__fastp(trimming_by_seqtype.ion.mix(trimming_illumina.bacteria)).trimmed

        trimmed_by_trimmomatic.mix(trimmed_by_fastp).branch {
            with_data: hasEnoughFastqData(it[1])
            insufficient_number_of_reads: true
        }
        .set { trimmed_branched }
        step_3TX_class__kraken(trimmed_branched.with_data)
    emit:
        no_reads = rawreads_branched.no_reads
        trimmed_with_data = trimmed_branched.with_data
        insufficient_number_of_reads = trimmed_branched.insufficient_number_of_reads
}

workflow {
    module_reads_processing(getInput())
}
