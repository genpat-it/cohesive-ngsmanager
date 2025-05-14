nextflow.enable.dsl=2

include { extractDsRef;getEmpty;flattenPath; parseMetadataFromFileName; executionMetadata; extractKey;taskMemory } from '../functions/common.nf'
include { isRunningFromSampleSheet } from '../functions/samplesheet.nf'
include { getSingleInput;getReferences;getReferenceCodes;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent;isSegmentedMapping } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '2AS_mapping'
def METHOD = 'minimap2' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def getExDt(reference, ex) {
    try {        
      reflist = getReferenceCodes()
      if (reflist.size() < 2) {
        return ex.dt
      } 
      if (!reference) {
        log.warn "${reference} should not be empty"
        return ex.dt
      }      
      if (!reflist.contains(reference)) {
        log.warn "${reference} not found in ${reflist}"
        return ex.dt
      }
      return ex.dt + reflist.indexOf(reference)
    } catch(Throwable t) {
        log.warn "getExDt: unexpected exception: ${t.asString()}"
        return ex.dt
    }
}

process minimap2 {
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_1"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    // memory { taskMemory( 10.GB, task.attempt ) }
    cpus 16    
    when:
      referencePath && referencePath.exists() && !referencePath.empty() && (isCompatibleWithSeqType(reads, ['nanopore'], task.process))
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(referencePath)
    output:
      path '*'
      tuple path("${base_ref}.sam"), val(reference), path(referencePath), emit: sam
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD,  [reference:reference, seq_type: 'nanopore'])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(reads.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      base_ref = "${base}_${reference.replace("_", "")}"      
      """
      minimap2 -ax map-ont ${referencePath} ${reads} -t ${task.cpus} -o ${base_ref}.sam
      """
}

process samtools {
    container "${LOCAL_REGISTRY}/bioinfo/samtools:0.1.19--f3869562fe"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    //memory { taskMemory( 6.GB, task.attempt ) }
      input:
      tuple path(sam), val(reference), path(referencePath)
    output:
      tuple path("${base_ref}_sorted.bam"), val(reference), emit: bam
      tuple path("${base_ref}.pileup"), val(reference), emit: pileup  
      path '*'
      path '*.sh', hidden: true
  
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*_sorted.bam*,*.pileup}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools.cfg" }
    script:
      md = parseMetadataFromFileName(sam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_minimap2_${reference}"
      """
      samtools view -bS -o ${base_ref}.bam ${sam} 2>> ${base_ref}.log
      samtools sort -@ 8 ${base_ref}.bam ${base_ref}_sorted 2>> ${base_ref}.log
      samtools index ${base_ref}_sorted.bam 2>> ${base_ref}.log
      samtools mpileup -d 1000 -A -Q 0 ${base_ref}_sorted.bam > ${base_ref}.pileup
	    """
}

process ivar {
    container "quay.io/biocontainers/ivar:1.3--h089eab3_1"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      tuple path(pileup), val(reference)
    output:
      path '*'
      path '*.sh', hidden: true
      tuple val(riscd), val(reference), path("${base_ref}.fasta"), emit: consensus    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_${METHOD}.cfg" }
    script:
      md = parseMetadataFromFileName(pileup.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${METHOD}_${reference}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
        # Note : samtools mpileup output must be piped into ivar consensus
        #cat ${pileup} | ivar consensus -p ${base_ref} -q 14 > ${base_ref}.log
        cat ${pileup} | ivar consensus -p ${base_ref} -q ${params.step_2AS_mapping__ivar___q} -m 1 > ${base_ref}.log
        mv ${base_ref}.fa ${base_ref}.fasta
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
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '{*.csv,*.check}'      
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_check.cfg" }
    script:
      //TODO fix output folder
      md = parseMetadataFromFileName(consensus.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
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
      base_ref = "${base}_minimap2_${reference}"
      """
      /scripts/coverage_plot.py ${coverage_depth} ${base_ref}_coverage_plot.png
	    """
}


workflow step_2AS_mapping__minimap2 {
  take: 
    reads
    reference 
  main:
    minimap2(reads, reference)
    pileup = samtools(minimap2.out.sam).pileup
    consensus = ivar(pileup).consensus

    coverage_minmax(samtools.out.bam, METHOD)
    coverage_minmax.out.coverage_depth | coverage_plot
    
    coverage = samtools_depth(samtools.out.bam, METHOD).coverage
    coverage.cross(consensus) { extractDsRef(it) }.map { 
        return [ it[0][0], it[1][2], it[0][1] ]
    }.set { coverageRefAndConsensus }
    coverageBasic = coverage_check(coverageRefAndConsensus, METHOD).coverage_basic

    crossedChecks = coverage_minmax.out.coverage_extra.cross(coverageBasic) { it[0] + "-" + it[1] }
    .map { [ it[0][0], it[0][1], it[0][2], it[1][2] ] }
    coverage_check_merge(crossedChecks, METHOD)

  emit:
    consensus = consensus
}


workflow {
    getSingleInput().cross(getReferences('any')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }
    step_2AS_mapping__minimap2(input.reads, input.refs)
}
