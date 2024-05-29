nextflow.enable.dsl=2

include { step_2AS_denovo__unicycler } from '../steps/step_2AS_denovo__unicycler'
include { step_4AN_AMR__blast } from '../steps/step_4AN_AMR__blast'
include { getSingleInput;getGenusSpeciesOptional } from '../functions/parameters.nf'
include { extractKey } from '../functions/common.nf'

workflow module_enterotoxin_saureus_finder {
    take: 
        trimmed
        genus_species
    main:
        assembly = step_2AS_denovo__unicycler(trimmed)

        assembly.cross(genus_species) { extractKey(it) }
            .multiMap { 
                assembly: it[0]
                species: it[1]
            }.set { assemblyAndSpecies }
        step_4AN_AMR__blast(assemblyAndSpecies.assembly, assemblyAndSpecies.species)
}

workflow {
    module_enterotoxin_saureus_finder(getSingleInput(), getGenusSpeciesOptional())
}
