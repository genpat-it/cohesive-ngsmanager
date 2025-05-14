nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;isSpeciesSupported } from '../functions/common.nf'
include { getInput;param } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_flaA'
def METHOD = 'flaA' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def MLST_SCHEMA_NAME = 'flaA'

def GENUS_ALLOWED = [
  'campylobacter'
]

process mlst_flaa {
    container "${LOCAL_REGISTRY}/bioinfo/mlst-w-db:2.23.0--60b8b2e3dd_231012.171302"
    containerOptions = "-v /bioinfonas:/bioinfonas:ro -v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro -u 0:0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 4.GB, task.attempt ) }
    cpus 8       
    when:
      isSpeciesSupported(genus_species, GENUS_ALLOWED, scaffolds200, task.process)
    input:
      tuple val(riscd_input), path(scaffolds200)
      val genus_species
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.csv,*.tsv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_mlst_flaa.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
        mlst --threads ${task.cpus} --scheme ${MLST_SCHEMA_NAME} ${scaffolds200} > ${base}_mbn_flaA.tsv 2> ${base}_mbn_flaA.log
        /scripts/process-mlst-result.py ${base}_mbn_flaA.tsv ${base}_mbn_flaa_cc.csv /NGStools/mlst/db/pubmlst
      """
}

workflow step_4TY_flaA__flaA {
    take: 
      assembly
      genus_species
    main:
      mlst_flaa(assembly, genus_species)
}

workflow {
    step_4TY_flaA__flaA(getInput(), param('genus_species'))
}