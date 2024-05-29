nextflow.enable.dsl=2

include { step_2AS_denovo__shovill } from '../steps/step_2AS_denovo__shovill'
include { getSingleInput } from '../functions/parameters.nf'

workflow module_wgs_bacteria {
    take: 
      trimmedReads
    main:
        step_2AS_denovo__shovill(trimmedReads)
    emit:
        step_2AS_denovo__shovill.out
}

workflow {
    module_wgs_bacteria(getSingleInput())
}