nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common'
include { getInput;isCompatibleWithSeqType;isIlluminaPaired } from '../functions/parameters.nf'
include { stepInputs;parseRISCD } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '0SQ_rawreads'
def METHOD = 'fastq' // XXX will be adjusted according to input RISCD
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process fastqc {
    container 'biocontainers/fastqc:v0.11.5_cv4'
    memory { taskMemory( 1.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      isCompatibleWithSeqType(reads, ['illumina_paired','ion'], null)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md_input, [dt: md_input.dt], md_input.acc, md_input.met, [seq_type:seq_type])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/result", pattern: '{*.zip,*.html}', saveAs: { filename -> filename.replaceFirst("-DT\\d+_", "-${ex.dt}_") }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/meta", pattern: '*.log', saveAs: { filename -> filename.replaceFirst("-DT\\d+_", "-${ex.dt}_") }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_${ex.dt}_fastqc.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      md_input = parseRISCD(riscd_input)       
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      seq_type = isIlluminaPaired(reads) ? 'illumina_paired' : 'ion'
      """
      fastqc $reads &> "${base}_fastqc.log" 
      """
}

process nanoplot {
    container 'quay.io/biocontainers/nanoplot:1.41.3--pyhdfd78af_0'
    memory { taskMemory( 4.GB, task.attempt ) }
    cpus { [8, params.max_cpus as int].min() }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      isCompatibleWithSeqType(reads, ['nanopore'], null)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      tuple val(riscd_input), path('*.txt'), emit: stats
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md_input, [dt: md_input.dt], md_input.acc, md_input.met, [seq_type:seq_type])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/result", pattern: '{*.txt,*.html}', saveAs: { filename -> filename.replaceFirst("-DT\\d+_", "-${ex.dt}_") }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/meta", pattern: '.command.log', saveAs: { "${base}_${ex.dt}_nanoplot.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_${ex.dt}_nanoplot.cfg" }
    script:
      md = parseMetadataFromFileName(reads.getName())
      md_input = parseRISCD(riscd_input)       
      seq_type = 'nanopore'
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      NanoPlot -t ${task.cpus} --fastq ${reads} --tsv_stats --no_static -o . -p ${base}_
      """
}

process nanopore_reads_check {
    container "quay.io/biocontainers/biopython:1.78"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    memory { taskMemory( 1.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(stats)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/result", pattern: '{*_readsCheck.csv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_${ex.dt}_readscheck.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md_input.ds}-${md_input.dt}_${md_input.met}/qc/meta", pattern: '.command.log', saveAs: { "${base}_${ex.dt}_readscheck.log" }
    script:
      md = parseMetadataFromFileName(stats.getName())
      md_input = parseRISCD(riscd_input)       
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      /scripts/SampleReadsCheck_nanopore.py -n $base -s $stats
      """
}

workflow step_0SQ_rawreads__fastq {
    take: data
    main:
      fastqc(data)
      nanoplot(data).stats | nanopore_reads_check
}

workflow {
  step_0SQ_rawreads__fastq(getInput())
}