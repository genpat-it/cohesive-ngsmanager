nextflow.enable.dsl=2

include { step_2AS_filtering__seqio } from '../steps/step_2AS_filtering__seqio'
include { step_3TX_species__vdabricate } from '../steps/step_3TX_species__vdabricate'
include { extractKey } from '../functions/common.nf'
include { getInput;getReference;getDS } from '../functions/parameters.nf'

workflow module_scaffolds_filtering {    
    take: 
        assembled
        reference 
        abricatedatabase 
    main:
        assembled.cross(abricatedatabase) { extractKey(it) }
            .map { 
                [ it[0][0], it[0][1], it[1][1] ] //scaffolds, db
            }.set { scaffoldsAndDatabase }
        calls = step_3TX_species__vdabricate(scaffoldsAndDatabase)
        calls
            .cross(assembled) { extractKey(it) } // [ [riscd, calls], [riscd, assembly] ]
            .cross(reference) { extractKey(it) } // [ [ [riscd, calls], [riscd, assembly] ], [ key, riscd, refid refpath ] ]
            .multiMap { 
                calls: it[0][0]
                assembly: it[0][1]
                reference: it[1][1..3]
            }.set { filt }
        step_2AS_filtering__seqio(filt.calls, filt.assembly, filt.reference)
}

workflow {
    module_scaffolds_filtering(getInput(), getReference('fa'), Channel.of([ getDS(), 'viruses_TREF' ]))
}