nextflow.enable.dsl=2

include { pipeline_ngsmanager } from './pipelines/pipeline_ngsmanager'
include { logHeader } from './functions/common.nf'

log.info logHeader('NGSMANAGER')

def helpMessage() {
    log.info"""
    Usage:

    # single-sample

    nextflow run http://gtc-devsrv:3000/bioinfo/ngsmanager --cmp 2021.TE.8600.1.96 --riscd 210108-11234826-2AS_mapping-bowtie

    Mandatory arguments
      --cmp                         sample id
      --riscd                       analysis result to be used as input
      --sample_type                 one of: 'bacterium', 'virus', 'sarscov2'

    # multi-sample

    nextflow run http://gtc-devsrv:3000/bioinfo/ngsmanager --samplesheet samples.csv 
    
    Mandatory arguments      
      --samplesheet                 samplesheet containing samples metadata

    """.stripIndent()
}

workflow {
    if (!params.containsKey('cmp') && !params.containsKey('samplesheet')){
        helpMessage()
        exit 1
    }
    pipeline_ngsmanager()
}