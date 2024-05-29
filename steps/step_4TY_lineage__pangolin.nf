nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_lineage'
def METHOD = 'pangolin' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process pangolin {
    container "${LOCAL_REGISTRY}/bioinfo/pangolin:v4.3.1-v0.1.12-v0.3.17-v1.21"
    containerOptions = "--user root"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 3.GB, task.attempt ) }
    input:
      tuple val(riscd_input), path(consensus)
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'lineage_report.csv', saveAs: { "${base}_lineage_report.csv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(consensus.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_pangolin"
      """
          pangolin ${consensus}  > ${base}.log
      """
}

workflow step_4TY_lineage__pangolin {
    take: 
      consensus
    main:
      pangolin(consensus)
}

workflow {
    step_4TY_lineage__pangolin(getInput())
}