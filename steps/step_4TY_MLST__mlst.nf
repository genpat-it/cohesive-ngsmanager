nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;param } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_MLST'
def METHOD = 'mlst' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process mlst {
    container "ghcr.io/genpat-it/mlst-w-db:2.23.0--60b8b2e3dd_231219.124455"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 4.GB, task.attempt ) }
    cpus { [8, params.max_cpus as int].min() }
    input:
      tuple val(riscd_input), path(assembly)
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.csv,*.tsv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      excluded_schemas = param('step_4TY_MLST__mlst___excluded_schemas')
      """
        mlst --threads ${task.cpus}  ${assembly} --exclude '${excluded_schemas}' > ${base}.tsv 2> ${base}.log
        /scripts/process-mlst-result.py ${base}.tsv ${base}_cc.csv /NGStools/mlst/db/pubmlst
      """
}

workflow step_4TY_MLST__mlst {
    take: 
      assembly
    main:
      mlst(assembly)
}

workflow {
  step_4TY_MLST__mlst(getInput())
}