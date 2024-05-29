nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4AN_AMR'
def METHOD = 'abricate'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process abricate {
    container 'staphb/abricate:1.0.0'
    memory { taskMemory( 250.MB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(assembly)
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_abricate.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.csv,*.summary,*calls.txt}'    
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """               
        abricate --db vfdb --csv ${assembly} > abricate.csv  
        abricate ${assembly} -db vfdb &>> ${base}_abricate.log >> ${base}_abricate_calls.txt        
        abricate --summary ${base}_abricate_calls.txt > ${base}_abricate.summary
      """
}

workflow step_4AN_AMR__abricate {
    take: data
    main:
      abricate(data);
}

workflow {
    step_4AN_AMR__abricate(getInput())
}
