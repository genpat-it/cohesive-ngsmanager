nextflow.enable.dsl=2

include { step_2AS_mapping__bowtie } from '../steps/step_2AS_mapping__bowtie'
include { step_2AS_mapping__ivar } from '../steps/step_2AS_mapping__ivar'
include { step_4AN_genes__prokka } from '../steps/step_4AN_genes__prokka'
include { step_4TY_lineage__pangolin } from '../steps/step_4TY_lineage__pangolin'
include { extractKey; getEmpty } from '../functions/common.nf'
include { getSingleInput;getReferenceOptional;getReference } from '../functions/parameters.nf'

def PROKKA_KINGDOM = 'Viruses'

workflow module_draft_genome {
    take: 
        reads
        reference 
        referenceGB
    main:
        reads.cross(reference) { extractKey(it) }.multiMap { 
            reads: it[0] // riscd, reads
            refs:  it[1][1..3] // riscd, code, path
        }.set { readsAndReferences }

        step_2AS_mapping__bowtie(readsAndReferences.reads, readsAndReferences.refs)

        consensus = step_2AS_mapping__ivar(readsAndReferences.reads, readsAndReferences.refs).consensus

        consensus.cross(referenceGB) { extractKey(it) }.map { 
            [ it[0][0], it[0][1], PROKKA_KINGDOM, it[1][1], it[1][2], it[1][3] ] // riscd assembly kingdom riscd_ref refid refpath]

        }.set { consensusKingdomReference }

        step_4AN_genes__prokka(consensusKingdomReference)
    }

workflow  {
    module_draft_genome(getSingleInput(), getReference('fa'), getReferenceOptional('gb'))
}