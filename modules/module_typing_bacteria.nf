nextflow.enable.dsl=2

include { step_2AS_mapping__bowtie } from '../steps/step_2AS_mapping__bowtie'
include { step_3TX_species__kmerfinder;getBacterialReferencePath } from '../steps/step_3TX_species__kmerfinder'
include { step_4AN_genes__prokka } from '../steps/step_4AN_genes__prokka'
include { step_4AN_AMR__abricate } from '../steps/step_4AN_AMR__abricate'
include { step_4AN_AMR__staramr } from '../steps/step_4AN_AMR__staramr'
include { step_4TY_cgMLST__chewbbaca } from '../steps/step_4TY_cgMLST__chewbbaca'
include { step_4TY_MLST__mlst } from '../steps/step_4TY_MLST__mlst'
include { step_4TY_flaA__flaA } from '../steps/step_4TY_flaA__flaA'
include { csv2map; extractKey; getEmpty } from '../functions/common.nf'
include { getTrimmedReads;getAssembly } from '../functions/parameters.nf'

workflow module_typing_bacteria {
    take: 
      trimmed
      assembly
    main:
      assigned_species = step_3TX_species__kmerfinder(assembly).assigned_species
      
      if (!params.skip_bestref_mapping) {
        trimmed.cross(assigned_species) { extractKey(it) }.multiMap { 
          trimmed: it[0]
          species: it[1][1]
          referencePath: it[1][2]
        }.set { trimAndAndSpecies }
        step_2AS_mapping__bowtie(trimAndAndSpecies.trimmed, trimAndAndSpecies.referencePath)
      } 

      step_4AN_AMR__abricate(assembly)

      step_4AN_genes__prokka(assembly.map{ [ it[0], it[1], 'Bacteria', '-', '-', getEmpty() ] })

      assembly.cross(assigned_species) { extractKey(it) }.multiMap { 
        assembly: it[0]
        species: it[1][1]
      }.set { assemblyAndSpecies }

      step_4AN_AMR__staramr(assemblyAndSpecies.assembly, assemblyAndSpecies.species)
      step_4TY_MLST__mlst(assemblyAndSpecies.assembly)
      step_4TY_flaA__flaA(assemblyAndSpecies.assembly, assemblyAndSpecies.species)
      step_4TY_cgMLST__chewbbaca(assemblyAndSpecies.assembly, assemblyAndSpecies.species, '')
    emit:
        genus_species = assigned_species
}

workflow {
    module_typing_bacteria(getTrimmedReads(true), getAssembly())
}