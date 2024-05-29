nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;stepInputs;extractKey } from '../functions/common.nf'
include { getSingleInput;isCompatibleWithSeqType;getHostReference;getRisCd } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '1PP_hostdepl'
def METHOD = 'minimap2' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process minimap2 {
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_1"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    // memory { taskMemory( 10.GB, task.attempt ) }
    cpus { [16, params.max_cpus as int].min() }
    when:
      isCompatibleWithSeqType(reads, ['nanopore'], task.process)
    input:
      tuple val(riscd_input), path(reads)
      tuple val(_), val(host_code), path(host)
    output:
      path '*'
      tuple path("${base_ref}.sam"), val(host_code), emit: sam
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD,  [reference:host_code, seq_type: 'nanopore'])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(reads.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      base_ref = "${base}_${host_code.replace("_", "")}"      
      """
      minimap2 -ax map-ont ${host} ${reads} -t ${task.cpus} -o ${base_ref}.sam
      """
}

process samtools {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    // memory { taskMemory( 18.GB, task.attempt ) }
    cpus { [16, params.max_cpus as int].min() }
    input:
      tuple path(sam), val(host_code)
    output:
      tuple val(riscd), path('*.fastq.gz'), emit: depleted
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "*.gz"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}_samtools.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_samtools.cfg" }
    script:
      md = parseMetadataFromFileName(sam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      base_ref = "${base}_${host_code.replace("_", "")}"
      riscd = getRisCd(md, ex, STEP, METHOD)                 
      """
        trap "rm '${base}_um.bam'" EXIT
		    samtools view -bS ${sam} -@ ${task.cpus} -o ${base}.bam 
        samtools view -b -f 4 ${base}.bam  -@ ${task.cpus} -o ${base}_um.bam  2>> ${base_ref}_samtools.log 
        samtools sort -n ${base}_um.bam -o ${base}.sorted -@ ${task.cpus} 2>> ${base_ref}_samtools.log 
        samtools bam2fq ${base}.sorted  | gzip > ${base_ref}.fastq.gz 2>> ${base_ref}_samtools.log
	    """
}

workflow step_1PP_hostdepl__minimap2 {
    take: 
      reads
      host
    main:
      minimap2(reads, host).sam | samtools
    emit:
      samtools.out.depleted 
}

workflow {
    getSingleInput().cross(getHostReference()) { extractKey(it) }
      .multiMap { 
          reads: it[0]
          host:  it[1]
      }.set { input }
  
    step_1PP_hostdepl__minimap2(input.reads, input.host)
}