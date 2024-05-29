nextflow.enable.dsl=2

include { module_denovo } from '../modules/module_denovo'
include { module_scaffolds_filtering } from '../modules/module_scaffolds_filtering'
include { module_draft_genome } from '../modules/module_draft_genome'
include { extractKey } from '../functions/common.nf'
include { getSingleInput;getHost;getReference;getReferenceOptional;getDS } from '../functions/parameters.nf'

workflow module_vdraft {
    take: 
        reads
        host
        reference
        referenceGB
        abricateDatabase
    main:
        denovoOut = module_denovo(reads, host);

        denovoOut.assembled
            .cross(reference) { extractKey(it) }
            .cross(abricateDatabase) { extractKey(it) }.multiMap { 
                assembly: it[0][0][0..1]
                reference: it[0][1]
                abricateDatabase: it[1]
            }.set { cARA }
        module_scaffolds_filtering(cARA.assembly, cARA.reference, cARA.abricateDatabase)
        
        denovoOut.depleted
            .cross(reference) { extractKey(it) }
            .cross(referenceGB) { extractKey(it) }
            .multiMap {
                depleted: it[0][0][0..1]
                reference: it[0][1]
                referenceGB: it[1]
            }
            .set { cDR }
        module_draft_genome(cDR.depleted, cDR.reference, cDR.referenceGB)
    }

workflow {
    module_vdraft(getSingleInput(), getHost(), getReference('fa'), getReferenceOptional('gb'), Channel.of([ getDS(), 'viruses_TREF' ]))
}
