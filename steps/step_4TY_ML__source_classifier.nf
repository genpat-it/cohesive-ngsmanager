nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_ML'
def METHOD = 'source_classifier'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def ALLELIC_PROFILE_ENCODING = 'crc32'

process source_classifier {
    container "nexus-prod.izs.intra:9091/bioinfo/ml_listeria_meat_classifier:1.0.0--907bd8d175_c19fb44f13"
    memory { taskMemory( 4.GB, task.attempt ) }
    cpus 16
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(profile)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true 
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'prediction.tsv', saveAs: { "${base}.tsv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(profile.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      """
      python /app/ml_listeria_meat_classifier.py --input ${profile} --model /app/xgb_15_0.2_ver250124.pkl --folder-output .
      """      
}

workflow step_4TY_ML__source_classifier {
    take: 
      allelic_profiles
    main:
      source_classifier(allelic_profiles)
}

workflow {
    step_4TY_ML__source_classifier(getInput())
}