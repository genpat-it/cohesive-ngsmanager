nextflow.enable.dsl=2

include { step_2AS_mapping__bowtie } from '../steps/step_2AS_mapping__bowtie'
include { step_3TX_species__kmerfinder;getBacterialReferencePath } from '../steps/step_3TX_species__kmerfinder'
include { step_4AN_genes__prokka } from '../steps/step_4AN_genes__prokka'
include { step_4AN_AMR__abricate } from '../steps/step_4AN_AMR__abricate'
include { step_4AN_AMR__staramr } from '../steps/step_4AN_AMR__staramr'
include { step_4TY_cgMLST__chewbbaca } from '../steps/step_4TY_cgMLST__chewbbaca'
include { step_4TY_wgMLST__chewbbaca } from '../steps/step_4TY_wgMLST__chewbbaca'
include { step_4TY_MLST__mlst } from '../steps/step_4TY_MLST__mlst'
include { step_4TY_flaA__flaA } from '../steps/step_4TY_flaA__flaA'
include { step_4TY_serotype__seqsero } from '../steps/step_4TY_serotype__seqsero'
include { module_brucella } from '../modules/module_brucella'
include { csv2map; extractKey; getEmpty } from '../functions/common.nf'
include { getTrimmedReads;getAssembly;optionalBoolean } from '../functions/parameters.nf'

workflow module_typing_bacteria {
    take: 
      trimmed
      assembly
    main:
      assigned_species = step_3TX_species__kmerfinder(assembly).assigned_species
      
      trimmed
        .cross(assigned_species) { extractKey(it) }
        .filter { !(it[1][1] ==~ /^Brucella_.*/) }
        .multiMap { 
          trimmed: it[0]
          species: it[1][1]
          referencePath: it[1][2]
        }.set { trimCrossSpecies }

      if (!params.skip_bestref_mapping) {
        step_2AS_mapping__bowtie(trimCrossSpecies.trimmed, trimCrossSpecies.referencePath)
      } 
      // specific processing for Brucella
      trimmed
        .cross(assigned_species) { extractKey(it) }
        .filter { it[1][1] ==~ /^Brucella_.*/ }
        .map { it[0] }.set { brucellas }
      module_brucella(brucellas)

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
      step_4TY_wgMLST__chewbbaca(assemblyAndSpecies.assembly, assemblyAndSpecies.species, '')
      step_4TY_serotype__seqsero(trimCrossSpecies.trimmed, trimCrossSpecies.species)
    emit:
        genus_species = assigned_species
}

workflow module_typing_bacteria_fasta {
    take: 
      assembly
    main:
      assigned_species = step_3TX_species__kmerfinder(assembly).assigned_species

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
      step_4TY_wgMLST__chewbbaca(assemblyAndSpecies.assembly, assemblyAndSpecies.species, '')
      // use fasta 
      step_4TY_serotype__seqsero(assemblyAndSpecies.assembly, assemblyAndSpecies.species)
}

workflow {
  if (optionalBoolean('module_typing_bacteria__use_fastq')) {
    module_typing_bacteria(getTrimmedReads(false), getAssembly())
  } else {
    module_typing_bacteria_fasta(getAssembly())
  }
}