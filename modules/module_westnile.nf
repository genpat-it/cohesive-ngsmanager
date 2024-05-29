nextflow.enable.dsl=2

// include { module_denovo } from '../modules/module_denovo'
// include { module_scaffolds_filtering } from '../modules/module_scaffolds_filtering'
include { step_4TY_lineage__westnile;getReferenceForLineage } from '../steps/step_4TY_lineage__westnile'
include { step_2AS_mapping__ivar } from '../steps/step_2AS_mapping__ivar'

include { extractKey } from '../functions/common.nf'
include { getSingleInput } from '../functions/parameters.nf'

workflow module_westnile {
    take: 
        reads        
    main:
        lineage = step_4TY_lineage__westnile(reads)

        reads.cross(lineage) { extractKey(it) }.multiMap { 
            reads: it[0]
            reference: getReferenceForLineage(it[1][1])
        }.set { readsAndRef }
        
        step_2AS_mapping__ivar(readsAndRef.reads, readsAndRef.reference).consensus
}

workflow {
    module_westnile(getSingleInput())
}
