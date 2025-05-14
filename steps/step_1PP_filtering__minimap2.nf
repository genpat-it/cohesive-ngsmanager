nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;stepInputs;extractKey } from '../functions/common.nf'
include { getSingleInput;isCompatibleWithSeqType;getReference;getRisCd} from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '1PP_filtering'
def METHOD = 'minimap2' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process minimap2 {
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_1"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    cpus { [16, params.max_cpus as int].min() }
    when:
      isCompatibleWithSeqType(reads, ['nanopore'], task.process)
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(referencePath)
    output:
      path '*'
      tuple path("${base}.sam"), val(reference), path(referencePath), emit: sam
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_ref], md, ex, STEP, METHOD, [reference:reference, seq_type: 'nanopore'])}' > ${base_ref}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base_ref}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      md = parseMetadataFromFileName(reads.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      base_ref = "${base}_${reference}"       
      """
      minimap2 -ax map-ont ${referencePath} ${reads} -t ${task.cpus} -o ${base}.sam
      """
}

process samtools {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 18.GB, task.attempt ) }
    cpus { [16, params.max_cpus as int].min() }
    input:
      tuple path(sam), val(reference), path(referencePath)
    output:
      tuple val(riscd), path('*.fastq.gz'), emit: filtered
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "*.gz"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*_sorted.bam*,*.vcf}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}_samtools.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_samtools.cfg" }
    script:
      md = parseMetadataFromFileName(sam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      base_ref = "${base}_${reference}"
      tmpName = "${base}_tmp"
      riscd = getRisCd(md, ex, STEP, METHOD)   
      """
		    samtools view -bS ${sam} -@ ${task.cpus} -o ${tmpName}.bam 
        samtools view -b -F 4 ${tmpName}.bam  -@ ${task.cpus} -o ${tmpName}_map.bam  2>> ${base_ref}_samtools.log
        samtools sort ${tmpName}_map.bam -o ${tmpName}.sorted.bam -@ ${task.cpus} 2>> ${base_ref}_samtools.log 
        samtools bam2fq ${tmpName}.sorted.bam | gzip > ${base_ref}.fastq.gz 2>> ${base_ref}_samtools.log
        samtools index ${tmpName}.sorted.bam 2>> ${base_ref}_samtools.log
	    """
}

workflow step_1PP_filtering__minimap2 {
    take: 
      reads
      reference
    main:
      minimap2(reads, reference).sam | samtools
    emit:
      samtools.out.filtered  
}

workflow {
    getSingleInput().cross(getReference('fa')) { extractKey(it) }
      .multiMap { 
          reads: it[0]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }
    step_1PP_filtering__minimap2(input.reads, input.refs)
}