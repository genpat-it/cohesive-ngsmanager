nextflow.enable.dsl=2

include { step_1PP_hostdepl__bowtie } from '../steps/step_1PP_hostdepl__bowtie'
include { step_2AS_mapping__bowtie } from '../steps/step_2AS_mapping__bowtie'
include { extractKey;getEmpty } from '../functions/common.nf'
include { getSingleInput;getHostOptional;getReference } from '../functions/parameters.nf'

workflow module_vdraft_light {
    take: 
        trimmedReads
        host
        reference
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
            .set { trimmedOrDepleted }

        trimmedOrDepleted.cross(reference) { extractKey(it) }
        .multiMap { 
            reads: it[0] // riscd, reads
            refs:  it[1][1..3] // riscd, code, path
        }.set { readsAndReferences }

        step_2AS_mapping__bowtie(readsAndReferences.reads, readsAndReferences.refs)
    }

workflow {
    module_vdraft_light(getSingleInput(), getHostOptional(), getReference('fa'))
}