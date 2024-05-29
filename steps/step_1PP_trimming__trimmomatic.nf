nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;extractKey } from '../functions/common'
include { getInput;isIonTorrent;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '1PP_trimming'
def METHOD = 'trimmomatic' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process trimmomatic {
    container 'quay.io/biocontainers/trimmomatic:0.36--6'
    memory { taskMemory( 1.GB, task.attempt ) }
        tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      tuple val(riscd), path('*.fastq.gz'), emit: fastq
      path '{*trimmomatic.log,*.json}'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fastq.gz'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*trimmomatic.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_trimmomatic.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      riscd = getRisCd(md, ex, STEP, METHOD)    
      if (isIlluminaPaired(reads)) { 
        """
          trimmomatic PE -threads 2 -phred33 $r1 $r2 ${base}_R1_trimmomatic.fastq.gz ${base}_R1_unpaired.fastq.gz ${base}_R2_trimmomatic.fastq.gz ${base}_R2_unpaired.fastq.gz ILLUMINACLIP:/usr/local/share/trimmomatic-0.36-6/adapters/NexteraPE-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:20:25 MINLEN:36 2>> ${base}_trimmomatic.log;
          cat ${base}_R1_unpaired.fastq.gz ${base}_R2_unpaired.fastq.gz > ${base}_unpaired_trimmomatic.fastq.gz
          rm  ${base}_R1_unpaired.fastq.gz ${base}_R2_unpaired.fastq.gz
        """
      } else if (isIonTorrent(reads)) {
        """
          trimmomatic SE -threads 2 -phred33 $r1 ${base}_R1_trimmomatic.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:55 2>> ${base}_trimmomatic.log
        """        
      }        
}

process sample_reads_check {
    container "quay.io/biocontainers/biopython:1.78"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    memory { taskMemory( 32.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(_), path(reads)
      tuple val(_), path(trimmed)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*SRC_raw.*[kg]'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '*SRC_raw.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*SRC_treads.*[kg]'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '{*SRC_treads.csv,*_readsCheck.csv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_readscheck.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      (t1,t2,u) = (trimmed instanceof java.util.Collection) ? trimmed : [trimmed, null, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      if (isIlluminaPaired(reads)) { 
        """
        /scripts/SampleReadsCheck.py -n $base -R1 $r1 -R2 $r2 -T1 $t1 -T2 $t2 -U $u > ${base}_SRC_raw.log
        cat ${base}_SRC_raw.log >	${base}_SRC_treads.log;
        """
      } else if (isIonTorrent(reads)) {
        """
        /scripts/SampleReadsCheck_ionTorrent.py -n $base -R1 $r1 -T1 $t1 > ${base}_SRC_raw.log
        """
      }
}

process fastqc {
    container 'biocontainers/fastqc:v0.11.5_cv4'
    memory { taskMemory( 400.MB, task.attempt ) }
        tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '{*.zip,*.html}', saveAs: { filename -> filename.replaceFirst("-DT\\d+_", "-${ex.dt}_") }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_fastqc.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      fastqc $reads > "${base}_fastqc.log" 2>&1
      """
}

workflow step_1PP_trimming__trimmomatic {
    take: rawreads
    main:
      trimmed = trimmomatic(rawreads).fastq;
      fastqc(trimmed)
      readsCheckInput = rawreads.cross(trimmed) { extractKey(it) }.multiMap { 
        rawreads: it[0]
        trimmed: it[1]
      }      
      sample_reads_check(readsCheckInput.rawreads, readsCheckInput.trimmed)
    emit:
      trimmed
}

workflow {
    step_1PP_trimming__trimmomatic(getInput())
}