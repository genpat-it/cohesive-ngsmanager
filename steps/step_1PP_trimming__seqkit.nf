nextflow.enable.dsl=2

include { parseMetadataFromFileName;getRisCd;executionMetadata;taskMemory;taskTime;extractKey;stepInputs } from '../functions/common.nf'
include { getInput;isCompatibleWithSeqType;param;optional } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '1PP_trimming'
def METHOD = 'seqkit' 

process seqkit {
    container "quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 3.GB, task.attempt ) }
    time { taskTime( 10.m, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, ['illumina_paired'], task.process)    
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      tuple val(riscd), path('*.fastq.gz'), emit: trimmed
      path '{*.sh,*.log}', hidden: true
       afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [adapterseq:adapter])}' > ${base}_input.json"
       publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.gz'
       publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.json,*.html}'
       publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
       publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      riscd = getRisCd(md, ex, STEP, METHOD)
      adapter = optional('step_1PP_trimming__seqkit__adapter')
        """
            seqkit grep -s -v -p ${adapter} ${r1} > ${base}_R1.temp.fastq.gz
            seqkit grep -s -v -p ${adapter} ${r2} > ${base}_R2.temp.fastq.gz
            seqkit pair -1 ${base}_R1.temp.fastq.gz -2 ${base}_R2.temp.fastq.gz
            mv ${base}_R1.temp.paired.fastq.gz ${base}_${adapter}_R1.fastq.gz
            mv ${base}_R2.temp.paired.fastq.gz ${base}_${adapter}_R2.fastq.gz
            rm -f ${base}_R1.temp.fastq.gz ${base}_R2.temp.fastq.gz
        """
}

workflow step_1PP_trimming__seqkit {
    take: 
      rawreads
    main:
      trimmed = seqkit(rawreads).trimmed
      readsCheckInput = rawreads.cross(trimmed) { extractKey(it) }.multiMap { 
        rawreads: it[0]
        trimmed: it[1]
      }      
    emit:
      trimmed  
}

workflow {
    step_1PP_trimming__seqkit(getInput())
}

