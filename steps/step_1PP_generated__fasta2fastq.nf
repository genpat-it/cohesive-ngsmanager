nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;flattenPath } from '../functions/common.nf'
include { getSingleInput } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '1PP_generated'
def METHOD = 'fasta2fastq' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process fasta2fastq {
    container "quay.io/biocontainers/biopython:1.78"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(assembly)
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'results/*.gz', saveAs: { filename -> flattenPath(filename) }      
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_fasta2fastq"
      """
        (mkdir input && mkdir results && cd input && ln -s ../${assembly} ${base}.fasta) 
        /scripts/fasta2fastq.py --fasta input --fastq results       
        gzip results/*.fastq
      """
}

workflow step_1PP_generated__fasta2fastq {
    take: 
      reads
    main:
      fasta2fastq(reads)
}

workflow {
    step_1PP_generated__fasta2fastq(getSingleInput())
}