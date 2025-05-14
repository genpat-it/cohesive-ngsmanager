nextflow.enable.dsl=2

include { step_2AS_mapping__snippy } from '../steps/step_2AS_mapping__snippy'
include { getSingleInput } from '../functions/parameters.nf'

def referenceCode = 'GCF_000740415.1'
def referencePath = "${params.assets_dir}/module_brucella/GCF_000740415.fasta"
def referenceRiscd = '220321-020220321135852203-2AS_import-external'

workflow module_brucella {
    take: 
        trimmed
    main:
        trimmed.multiMap {
            trimmed: it
            reference: [ referenceRiscd, referenceCode, file(referencePath) ]
        }.set { trAndRef }
        step_2AS_mapping__snippy(trAndRef.trimmed, trAndRef.reference)
    }

workflow {
    module_brucella(getSingleInput())
}