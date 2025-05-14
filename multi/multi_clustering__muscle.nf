nextflow.enable.dsl=2

include { taskMemory } from '../functions/common.nf'
include { optionalBoolean;param;optionalOrDefault;getInput } from '../functions/parameters.nf'
include { muscle;prepare_multifasta;constant_sites;seqkit_stats;amas_summary } from '../components/alignment.nf'
include { multi_alignment_analysis } from '../components/subworkflows.nf'

workflow multi_clustering__muscle {
    take: 
        input
        raw_metadata
        geodata
        nomenclature
    main:
        summary = input.collectFile { [ "summary.tsv", it[0] + '\t' + it[1][0].getName() + '\n' ] }
        files= input.collect { it[1] }

        multifasta = prepare_multifasta(summary, files).multifasta
        
        seqkit_stats(multifasta)
        ma = muscle(multifasta).alignment
        amas_summary(ma)
     
        fconst = constant_sites(ma).constant_sites

        multi_alignment_analysis(ma, raw_metadata, geodata, nomenclature, false, fconst)
}

workflow {
  multi_clustering__muscle(getInput(),  param('metadata'), param('geodata'), optionalOrDefault('reportree__nomenclature', '/tmp/fake'));
}