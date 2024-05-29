nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;taskTime } from '../functions/common.nf'
include { getInput;isCompatibleWithSeqType;paramWrap;optWrap } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '1PP_trimming'
def METHOD = 'chopper' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process chopper {
    container "quay.io/biocontainers/chopper:0.7.0--hdcf5f25_0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 5.GB, task.attempt ) }
    time { taskTime( 10.m, task.attempt ) }
    cpus { [16, params.max_cpus as int].min() }       
    when:
      isCompatibleWithSeqType(reads, ['nanopore'], task.process)    
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      tuple val(riscd), path('*.fastq.gz'), emit: trimmed
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [seq_type: 'nanopore'])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.gz'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.json,*.html}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      quality = paramWrap('step_1PP_trimming__chopper__quality', '--quality {}')
      minlength = paramWrap('step_1PP_trimming__chopper__minlength', '--minlength {}')
      maxlength = optWrap('step_1PP_trimming__chopper__maxlength', '--maxlength {}')
      headcrop = optWrap('step_1PP_trimming__chopper__headcrop', '--headcrop {}')
      tailcrop = optWrap('step_1PP_trimming__chopper__tailcrop', '--tailcrop {}')
        """
          zcat $r1 | chopper ${quality} ${minlength} ${maxlength} \
          ${headcrop} ${tailcrop} \
          --threads $task.cpus | gzip > ${base}_R1.fastq.gz
        """              
}

process nanoplot {
    container 'quay.io/biocontainers/nanoplot:1.41.3--pyhdfd78af_0'
    memory { taskMemory( 4.GB, task.attempt ) }
    cpus { [8, params.max_cpus as int].min() }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '{*.txt,*.html}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.log', saveAs: { "${base}_nanoplot.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_nanoplot.cfg" }
    script:
      md = parseMetadataFromFileName(reads.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      NanoPlot -t ${task.cpus} --fastq ${reads} --tsv_stats --no_static -o . -p ${base}_
      """
}

workflow step_1PP_trimming__chopper {
    take: 
      rawreads
    main:
      trimmed = chopper(rawreads).trimmed
      nanoplot(trimmed)
    emit:
      trimmed      
}

workflow {
    step_1PP_trimming__chopper(getInput())
}