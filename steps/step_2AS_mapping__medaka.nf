nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;extractKey;taskMemory;stepInputs;getRisCd;extractDsRef } from '../functions/common.nf'
include { getSingleInput;getReference;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '2AS_mapping'
def METHOD = 'medaka' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process medaka {
    container "quay.io/biocontainers/medaka:1.8.0--py38hdaa7744_0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    // memory { taskMemory( 8.GB, task.attempt ) }
    cpus { [32, params.max_cpus as int].min() }
    when:
      reference_path && reference_path.exists() && !reference_path.empty() && (isCompatibleWithSeqType(reads, ['nanopore'], task.process))
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(reference_path)
    output:
      path '*_input.json'
      tuple val(riscd), path("${base_ref}.fasta"), emit: consensus    
      tuple path("${base_ref}.bam"), val(reference), emit: bam    
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_ref], md, ex, STEP, METHOD, [reference:reference])}' > ${base_ref}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "*.fasta"    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: "*.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base_ref}.log" }    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      md = parseMetadataFromFileName(reads.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${METHOD}_${reference}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
        #trap "find -name "*.sam" -delete" EXIT
        cp ${reference_path} _${reference_path}
        medaka_consensus \
            -t ${task.cpus} \
            -i ${reads} \
            -d _${reference_path} \
            -o ./
        mv consensus.fasta ${base_ref}.fasta
        mv calls_to_draft.bam ${base_ref}.bam        
      """      
}


process coverage_minmax {
    container "${LOCAL_REGISTRY}/bioinfo/samtools:0.1.19--f3869562fe"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 4.GB, task.attempt ) }
        input:
      tuple path(bam), val(reference)
      val(method)
    output:
      path '*.csv'
      path '*.sh', hidden: true
      tuple path("${base_ref}_samtools_depth.txt"), val(reference), emit: coverage_depth  
      tuple val(md.ds), val(reference), path("${base_ref}_import_coverage_minmax.csv"), emit: coverage_extra
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_minmax.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${method}_${reference}"
      """
      samtools view -F 4 -c ${bam} > samtools_view.txt
      samtools depth ${bam} > ${base_ref}_samtools_depth.txt
      /scripts/coverage_minmax.py ${md.cmp} ${md.ds} samtools_view.txt ${base_ref}_samtools_depth.txt ${base_ref}_import_coverage_minmax.csv
	    """
}

process samtools_depth {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
        input:
      tuple path(bam), val(reference)
      val(method)
    output:
      tuple path("${base_ref}.coverage"), val(reference), emit: coverage
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools_depth.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${method}_${reference}"
      """
      samtools depth -a ${bam} | awk '{ if (\$3!=0) c++;s+=\$3}{h++} END { if (c!=0) print s/c; else print 0;if (h!=0) print c/h; else print 0 }' > ${base_ref}.coverage
	    """
}

process coverage_check {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      tuple path(coverage), path(consensus), val(reference)
      val(context)
    output:
      tuple val(md.ds), val(reference), path("${base_ref}_import_coverage.csv"), emit: coverage_basic
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.csv,*.check}'      
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_check.cfg" }
    script:
      //TODO fix output folder
      md = parseMetadataFromFileName(consensus.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${context}_${reference}"
      coverages= coverage.toRealPath().toFile().readLines()
      """
      /scripts/coverage.py ${coverage} ${consensus} ${reference} ${md.cmp} ${md.ds.substring(2)} ${base_ref}.check ${base_ref}_import_coverage.csv 
	    """
}

process coverage_check_merge {
  container "ubuntu:20.04"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory { taskMemory( 200.MB, task.attempt ) }
  input:
    tuple val(key), val(reference), path(covMinMax), path(covBasic)
    val(method)
  output:
    tuple val(md.ds), path("${base_ref}_import_coverage_merged.csv"), emit: coverage_merged
    path '*'
    path '*.sh', hidden: true
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.csv'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_check_merge.cfg" }
  script:
    md = parseMetadataFromFileName(covMinMax.getName())
    base = "${md.ds}-${ex.dt}_${md.cmp}"
    base_ref = "${base}_${method}_${reference}"
    """
    paste -d, ${covMinMax} ${covBasic} | cut -d, -f1,2,3,4,5,8,9,10,11,12,13 > ${base_ref}_import_coverage_merged.csv
    """  
}

process coverage_plot {
    container "quay.io/biocontainers/matplotlib:3.1.2"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    input:
      tuple path(coverage_depth), val(reference)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.png'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_plot.cfg" }
    script:
      md = parseMetadataFromFileName(coverage_depth.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_bowtie_${reference}"
      """
      /scripts/coverage_plot.py ${coverage_depth} ${base_ref}_coverage_plot.png
	    """
}

workflow step_2AS_mapping__medaka {
  take: 
    reads
    reference 
  main:
    medaka(reads, reference)

    coverage_minmax(medaka.out.bam, METHOD)
    coverage_minmax.out.coverage_depth | coverage_plot
    
    coverage = samtools_depth(medaka.out.bam, METHOD).coverage
    coverage.cross(medaka.out.consensus) { extractDsRef(it) }.map { 
        return [ it[0][0], it[1][1], it[0][1] ]
    }.set { coverageRefAndConsensus }
    coverageBasic = coverage_check(coverageRefAndConsensus, METHOD).coverage_basic

    crossedChecks = coverage_minmax.out.coverage_extra.cross(coverageBasic) { it[0] + "-" + it[1] }
    .map { [ it[0][0], it[0][1], it[0][2], it[1][2] ] }
    coverage_check_merge(crossedChecks, METHOD)

  emit:
    consensus = medaka.out.consensus
}

workflow {
    getSingleInput().cross(getReference('fa')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }
    step_2AS_mapping__medaka(input.reads, input.refs)
}