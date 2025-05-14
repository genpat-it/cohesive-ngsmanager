nextflow.enable.dsl=2

include { step_0SQ_preprocess__fastq } from '../steps/step_0SQ_preprocess__fastq'
include { step_0SQ_preprocess__bcl2fastq } from '../steps/step_0SQ_preprocess__bcl2fastq'
include { getRawReadsFromSampleSheet } from '../functions/parameters.nf'

workflow module_reads_preprocessing {
    take: 
      rawReads
    main:
      preProcessedReads = step_0SQ_preprocess__fastq(rawReads)
    emit:
      preProcessedReads
}

workflow {
    module_reads_preprocessing(getRawReadsFromSampleSheet())
}