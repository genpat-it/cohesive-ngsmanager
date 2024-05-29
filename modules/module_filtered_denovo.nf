nextflow.enable.dsl=2

include { step_1PP_filtering__bowtie } from '../steps/step_1PP_filtering__bowtie'
include { step_2AS_denovo__spades } from '../steps/step_2AS_denovo__spades'
include { extractKey } from '../functions/common.nf'
include { getSingleInput;getReference } from '../functions/parameters.nf'

workflow module_filtered_denovo {
    take: 
        reads
        reference 
    main:
        reads.cross(reference) { extractKey(it) }.multiMap { 
            reads: it[0] // riscd, reads
            refs:  it[1][1..3] // riscd, code, path
        }.set { readsAndReferences }

        filtered = step_1PP_filtering__bowtie(readsAndReferences.reads, readsAndReferences.refs)
        assembled = step_2AS_denovo__spades(filtered)
    emit:
        assembled
        filtered
}

workflow {
    module_filtered_denovo(getSingleInput(), getReference('fa'))
}