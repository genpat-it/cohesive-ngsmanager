nextflow.enable.dsl=2

include { extractKey } from '../functions/common.nf'
include { step_2AS_mapping__ivar } from '../steps/step_2AS_mapping__ivar'
include { getSingleInput;getReferences } from '../functions/parameters.nf'

workflow module_segmented {
    take: 
        reads
        reference 
    main:
        step_2AS_mapping__ivar(reads, reference)
    emit:
        step_2AS_mapping__ivar.out.consensus
}

workflow {
    getSingleInput().cross(getReferences('any')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }

    module_segmented(input.reads, input.refs)
}