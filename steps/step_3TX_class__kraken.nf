nextflow.enable.dsl=2

include { getEmpty;flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { isSarsCov2;isPositiveControlSarsCov2;isNegativeControlSarsCov2;isNGSMG16S } from '../functions/samplesheet'
include { isFullOutput;getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters'
include { stepInputs;getRisCd } from '../functions/common.nf'

def db_kraken="/biowork/databases/PROGRAMS/kraken/minikraken_20171019_8GB/"
def db_bracken="/biowork/databases/PROGRAMS/Bracken-master/minikraken_8GB_125mers_distrib.txt"

def ex = executionMetadata()

def STEP = '3TX_class'
def METHOD = 'kraken' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process kraken {
    container "quay.io/biocontainers/kraken:1.0--pl5.22.0_0"
    containerOptions = "-v /biowork:/biowork:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 8.GB, task.attempt ) }
    when:
      (isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)) && !(isSarsCov2(riscd_input) || isPositiveControlSarsCov2(riscd_input) || isNegativeControlSarsCov2(riscd_input) || isNGSMG16S(riscd_input))
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path("${base}_kraken.tsv"), emit: report  
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_kraken.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.tsv,*.txt}'
    script:
      (t1,t2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      options = isIlluminaPaired(reads) ? ' --paired ' : ''
      options += t1.getName().endsWith(".gz") ? " --gzip-compressed " : ""
      """
        kraken --db ${db_kraken} --threads 8 --fastq-input ${options} ${t1} ${t2} > ${base}_kraken.txt 2> ${base}_kraken.log
        kraken-report --db ${db_kraken} ${base}_kraken.txt > ${base}_kraken.tsv
        ${isFullOutput()} || rm *.txt
      """    
}

process braken {
    container "quay.io/biocontainers/bracken:1.0.0--1"
    containerOptions = "-v /biowork:/biowork:ro -v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 3.GB, task.attempt ) }
    input:
      path(kraken_report)
    output:
      path '*'
      tuple val(riscd), path("${base}_bracken_genus.tsv"), emit: genus_report
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*_import_taxa.csv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_bracken.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.tsv'
    script:
      md = parseMetadataFromFileName(kraken_report.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
        est_abundance.py -i ${kraken_report} -k ${db_bracken} -o ${base}_bracken_genus.tsv -l G -t 10 > ${base}_bracken_genus.log
        est_abundance.py -i ${kraken_report} -k ${db_bracken} -o ${base}_bracken_species.tsv -l S -t 10 > ${base}_bracken_species.log
        /scripts/create_import.py -ds ${md.ds} -o ${base}_import_taxa.csv -g ${base}_bracken_genus.tsv -sp ${base}_bracken_species.tsv -k ${kraken_report}
      """    
}

workflow step_3TX_class__kraken {
    take: reads
    main:
      kraken(reads).report | braken
     emit:
       genus_report = braken.out.genus_report
}

workflow {
  step_3TX_class__kraken(getSingleInput())
}
