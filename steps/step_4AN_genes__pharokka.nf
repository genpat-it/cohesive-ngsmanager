nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;isSpeciesSupported } from '../functions/common.nf'
include { getInput;param } from '../functions/parameters.nf'
include { flattenPath;stepInputs;getRisCd } from '../functions/common.nf'

def PHAROKKA_DB=param('step_4AN_genes__pharokka___db')

def ex = executionMetadata()

def STEP = '4AN_genes'
def METHOD = 'pharokka' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process pharokka {
    container "quay.io/biocontainers/pharokka:1.7.2--pyhdfd78af_0"
    containerOptions = "-v ${PHAROKKA_DB}:${PHAROKKA_DB}:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    cpus 64
    // memory { taskMemory( 62.GB, task.attempt ) }
    input:
      tuple val(riscd_input), path(assembly)
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{output/*pharokka.gff,output/*pharokka.gbk,output/*final_merged_output.tsv}', saveAs: { filename -> flattenPath(filename) }           

    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      predictor = param('step_4AN_genes__pharokka___predictor')
      """
        pharokka.py -i ${assembly} -o output -p ${base} -d ${PHAROKKA_DB} -g ${predictor} -t ${task.cpus}
      """    
}

workflow step_4AN_genes__pharokka {
    take:
      assembly
    main:
      pharokka(assembly)
}

workflow {
  step_4AN_genes__pharokka(getInput())
}