nextflow.enable.dsl=2

include {  parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getGenusSpeciesOptional;getInput } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4AN_AMR'
def METHOD = 'blast' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def blastDBs = [
  staphylococcus_aureus : "/biowork/databases/PROGRAMS/BLAST/database_s.aureus_virulence.fsa"
]

process blast {
    container 'ncbi/blast:2.12.0'
    containerOptions = "-v /biowork:/biowork:ro"
    memory { taskMemory( 4.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      blastDBs.containsKey(genus_species?.toLowerCase())    
    input:
      tuple val(riscd_input), path(scaffolds200)
      tuple val(_), val(genus_species)
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.tsv'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_blast.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      blastDB = blastDBs.get(genus_species?.toLowerCase())
      """                
        blastn -query ${scaffolds200} -db ${blastDB} -outfmt 6 -out ${base}_blast_out.tsv      
      """
}

workflow step_4AN_AMR__blast {
    take: 
      assembly
      species
    main:
      blast(assembly, species);
}

workflow {
    step_4AN_AMR__blast(getInput(), getGenusSpeciesOptional())
}

