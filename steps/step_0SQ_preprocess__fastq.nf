nextflow.enable.dsl=2

include { taskMemory;executionMetadata } from '../functions/common.nf'
include { getDsMetadata } from '../functions/samplesheet.nf'
include { getRisCd } from '../functions/common.nf'

def extractDs(fileName) {
  def matcher = (fileName =~ /(DS\d+)_.+/)
  if (matcher.matches())  {
      return matcher.group(1)
  } else {
    log.error "Wrong filename pattern for raw reads: ${fileName}"
    return ""
  }
}

def ex = executionMetadata()

def STEP = '0SQ_rawreads'
def METHOD = 'fastq' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process import_fastq {
    container 'ubuntu:20.04'
    tag "${cmp}/${ds}/${ex?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    errorStrategy 'ignore'
    stageInMode 'symlink'
    maxForks 1
    input:
      tuple val(ds), path(r1_lanes), path(r2_lanes), path(samplesheet)
    output:
      tuple val(riscd), path('*.gz'), emit: fastq
      path '*.csv'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${anno}/${cmp}/${STEP}/${ds}-${ex.dt}_${METHOD}/result", pattern: "*.gz"
    publishDir mode: 'rellink', "${params.outdir}/${anno}/${cmp}/${STEP}/${ds}-${ex.dt}_${METHOD}/meta", pattern: '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${anno}/${cmp}/${STEP}/${ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_preprocessing.cfg" }
    script:
      sampleSheetMetaData = getDsMetadata(ds)
      cmp = sampleSheetMetaData.SAMPLE_PROJECT.replaceAll("-", ".")
      anno = cmp.substring(0,4)
      base = "${ds}-${ex.dt}_${cmp}"
      riscd = getRisCd([ds:ds], ex, STEP, METHOD)   
      if (!r2_lanes) {
        """
        mv ${r1_lanes} ${ds}-${ex.dt}_${cmp}_R1.fastq.gz
        grep "SAMPLE_ID," ${samplesheet} > ${base}_samplesheet.csv  && grep "^${cmp}," ${samplesheet} >> ${base}_samplesheet.csv 
        """
      } else if (r1_lanes instanceof java.util.Collection) {
        """
        cat ${r1_lanes} > ${ds}-${ex.dt}_${cmp}_R1_001.fastq.gz
        cat ${r2_lanes} > ${ds}-${ex.dt}_${cmp}_R2_001.fastq.gz
        grep "SAMPLE_ID," ${samplesheet} > ${base}_samplesheet.csv  && grep "^${cmp}," ${samplesheet} >> ${base}_samplesheet.csv 
        """
      } else {
        """
        mv ${r1_lanes} ${ds}-${ex.dt}_${cmp}_R1_001.fastq.gz
        mv ${r2_lanes} ${ds}-${ex.dt}_${cmp}_R2_001.fastq.gz
        grep "SAMPLE_ID," ${samplesheet} > ${base}_samplesheet.csv  && grep "^${cmp}," ${samplesheet} >> ${base}_samplesheet.csv 
        """
      }     
}

workflow step_0SQ_preprocess__fastq {
    take: data
    main:
      import_fastq(data)
    emit:
      import_fastq.out.fastq
}