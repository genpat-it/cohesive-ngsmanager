nextflow.enable.dsl=2

include {  flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

FILTERABLE_REFERENCES_PATH ='/databases/REFERENCES/virus'

def ex = executionMetadata()

def STEP = '2AS_denovo'
def METHOD = 'spades' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def isReferenceFilterable(refCode, _path) {
    try {
      def referencePath = (_path instanceof java.util.Collection) ? _path.flatten()[0] : _path
      if (referencePath.empty) {
        //it means that no reference is provided
        log.warn "no reference provided - skipping scaffold filtering"
        return false
      }
      println "${referencePath.toRealPath().toString()} vs ${FILTERABLE_REFERENCES_PATH}"
      //we can filter a reference only if present in abricate virus DB
      //assuming: reference path in  "referencesDir" -> then it is filterable
      def filterable = referencePath.toRealPath().toString().contains(FILTERABLE_REFERENCES_PATH)
      if (!filterable) {
        log.warn "a reference can be filtered only if present in our abricate virus DB - skipping scaffold filtering for '${refCode}'"
      }
      return filterable
    } catch(Throwable t) {
        log.error "error while executing 'isReferenceFilterable' on '${referencePath}' (${t.toString()})"
        return false
    } 
}

process filter {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 250.MB, task.attempt ) }
    when:
       isReferenceFilterable(reference, referencePath)
    input:
      tuple val(riscd_input), path(calls)
      tuple val(riscd_input2),  path(l200)
      tuple val(riscd_ref), val(reference), path(referencePath)
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_input2, riscd_ref], md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_filter.cfg" }
    script:
      md = parseMetadataFromFileName(calls.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      output = l200?.getName().replace(".fasta", "_${reference.replace('_', '')}.fasta")
      """
      /scripts/cleanDenovo.py ${l200} ${calls} ${reference} ${output}
      """
}

workflow step_2AS_filtering__seqio {
    take: 
      calls
      assembly
      reference
    main:
      filter(calls, assembly, reference)
}
