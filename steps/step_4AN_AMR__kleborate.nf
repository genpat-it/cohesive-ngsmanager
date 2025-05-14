nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory} from '../functions/common.nf'
include { getSingleInput;optWrap;optional;param;optionalBoolean } from '../functions/parameters.nf'
include { stepInputs;flattenPath } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4AN_AMR'
def METHOD = 'kleborate' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process kleborate {
    container "nexus-prod.izs.intra:9091/bioinfo/kleborate:2.1--2a44e8eeaf"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      tuple val(riscd_input), path(assembly)
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [method: optionalprofile])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "${base}/*.txt" , saveAs: { "${base}.txt" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script: 
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      optionalprofile = optional('step_4AN_AMR__kleborate__profile')
      """
        /usr/local/bin/_entrypoint.sh kleborate -a ${assembly} -o ${base} ${optionalprofile} --trim_headers
      """
}

workflow step_4AN_AMR__kleborate {
    take: 
      assembly
    main:
      kleborate(assembly)
}

workflow {
    step_4AN_AMR__kleborate(getSingleInput())
}