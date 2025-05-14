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

process import_bcl2fastq {
    container 'ubuntu:20.04'
    tag "${cmp}/${ds}/${ex?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    errorStrategy 'ignore'
    maxForks 1
    input:
      tuple path(r1), path(r2) ,path(samplesheet) //DS10178196-DT200310_2020.TE.4395.1.10_R1_001.fastq.gz
    output:
      tuple val(riscd), path('*.gz'), emit: fastq
      path '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${anno}/${cmp}/${STEP}/${ds}-${ex.dt}_${METHOD}/result", pattern: '*.gz'
    publishDir mode: 'rellink', "${params.outdir}/${anno}/${cmp}/${STEP}/${ds}-${ex.dt}_${METHOD}/meta", pattern: '*.csv'
    when:
      extractDs(r1.getName())
    script:
      ds = extractDs(r1.getName())
      sampleSheetMetaData = getDsMetadata(ds)
      cmp = sampleSheetMetaData.SAMPLE_PROJECT.replaceAll("-", ".")
      anno = cmp.substring(0,4)
      base = "${ds}-${ex.dt}_${cmp}"    
      r1Suffix = r1.getName() -~ /DS\d+/
      r2Suffix = r2.getName() -~ /DS\d+/
      riscd = getRisCd([ds:ds], ex, STEP, METHOD)      
      """
      mv ${r1} ${ds}-${ex.dt}_${cmp}${r1Suffix}
      mv ${r2} ${ds}-${ex.dt}_${cmp}${r2Suffix}
      head -n 2 ${samplesheet} > SampleSheetCmp.csv
      grep "^${cmp}," ${samplesheet} >> SampleSheetCmp.csv 
      """
}

workflow step_0SQ_preprocess__bcl2fastq {
    take: data
    main:
      import_bcl2fastq(data)
    emit:
      import_bcl2fastq.out.fastq
}