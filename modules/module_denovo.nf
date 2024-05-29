nextflow.enable.dsl=2

include { step_1PP_hostdepl__bowtie } from '../steps/step_1PP_hostdepl__bowtie'
include { step_2AS_denovo__spades } from '../steps/step_2AS_denovo__spades'
include { extractKey; getEmpty } from '../functions/common.nf'
include { getSingleInput;getHost } from '../functions/parameters.nf'

workflow module_denovo {
    take: 
        trimmedReads
        host 
    main:
        trimmedReads.cross(host) { extractKey(it) }
            .map { [ it[0][0], it[0][1], it[1][1] ] } //riscd, reads, host
            .branch {
                with_host: it[1][1]
                without_host: true
            }
        .set { branchedTrimmed }

        depleted = step_1PP_hostdepl__bowtie(branchedTrimmed.with_host)

        branchedTrimmed.without_host
            .mix(depleted)
            .map { it[0,1] }
            .set{ denovoInput }
        assembled = step_2AS_denovo__spades(denovoInput)
    emit:
        assembled
        depleted
}

workflow {
    module_denovo(getSingleInput(), getHost())
}