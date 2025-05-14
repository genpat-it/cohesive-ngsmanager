nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory} from '../functions/common.nf'
include { getSingleInput;optWrap;optional;param;optionalBoolean } from '../functions/parameters.nf'
include { stepInputs;flattenPath } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4AN_AMR'
def METHOD = 'amrfinderplus' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process amrfinderplus {
    container "ncbi/amr"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      tuple val(riscd_input), path(assembly)
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [plus: optionalplus, taxon:optionaltaxon])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.tsv', saveAs: { filename -> "${base}.tsv" }  
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script: 
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      optionaltaxon = optional('step_4AN_AMR__amrfinderplus__taxon')
      optionalplus = optionalBoolean('step_4AN_AMR__amrfinderplus__plus') ? "--plus" : ""
      if (optionaltaxon == '') {
        """
          amrfinder -n ${assembly} -o ${base}.tsv ${optionaltaxon} ${optionalplus}
        """
      } else {    
        """
          amrfinder -n ${assembly} -o ${base}.tsv --organism ${optionaltaxon} ${optionalplus}
        """
      }

}

workflow step_4AN_AMR__amrfinderplus {
    take: 
      assembly
    main:
      amrfinderplus(assembly)
}

workflow {
    step_4AN_AMR__amrfinderplus(getSingleInput())
}