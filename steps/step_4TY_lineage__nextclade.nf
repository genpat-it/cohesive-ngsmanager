nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;isSpeciesSupported } from '../functions/common.nf'
include { getSingleInput;param } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_lineage'
def METHOD = 'nextclade' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_SPECIES_ALLOWED = [
  "Virus_Dengue"
]

process nextclade {
  container "${LOCAL_REGISTRY}/bioinfo/nextclade:3.10.2--ac8af1a80a"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory { taskMemory( 8.GB, task.attempt ) }
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.tsv}'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.json,*.md}'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
  afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [dataset: "community/v-gen-lab/dengue/denv${serotype}"])}' > ${base}_input.json"
  input:
    tuple val(riscd_input), path(consensus)
    val(genus_species)
    val(serotype)
  output:
    path '*'
    path '{*.sh,*.log}', hidden: true
  when:
    isSpeciesSupported(genus_species, GENUS_SPECIES_ALLOWED, consensus, task.process) &&  (serotype instanceof Number || serotype.isNumber())
  script:
    md = parseMetadataFromFileName(consensus.getName())
    base = "${md.ds}-${ex.dt}_${md.cmp}_nextclade_denv${serotype}"
    """
      nextclade run --input-dataset /datasets/denv${serotype} --output-tsv=${base}.tsv ${consensus} 
      cp /datasets/denv${serotype}/CHANGELOG.md ${base}_changelog.md
    """
}

workflow step_4TY_lineage__nextclade {
    take: 
      consensus
      genus_species
      serotype
    main:
      nextclade(consensus, genus_species, serotype)
}

workflow {    
    step_4TY_lineage__nextclade(getSingleInput(), param('genus_species'), param('serotype'))
}