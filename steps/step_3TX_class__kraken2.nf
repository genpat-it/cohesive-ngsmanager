nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { isSarsCov2;isPositiveControlSarsCov2;isNegativeControlSarsCov2 } from '../functions/samplesheet'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters'
include { stepInputs;getRisCd } from '../functions/common.nf'

def KRAKEN2_DB="/bioinfonas/databases/PROGRAMS/kraken2/k2_pluspf_20210517"

def ex = executionMetadata()

def STEP = '3TX_class'
def METHOD = 'kraken2' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process kraken2 {
    container "quay.io/biocontainers/kraken2:2.1.2--pl5262h7d875b9_0"
    containerOptions = "-v /bioinfonas/databases:/bioinfonas/databases:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 62.GB, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path("${base}_kraken.tsv"), emit: report  
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_kraken.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.tsv,*.gz}'
    script:
      (t1,t2) = reads
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      options = isIlluminaPaired(reads) ? ' --paired ' : ''
      options += t1.getName().endsWith(".gz") ? " --gzip-compressed " : ""
      """
        kraken2 --db ${KRAKEN2_DB} --threads 64 ${options} ${t1} ${t2} --output ${base}_kraken.txt --report ${base}_kraken.tsv 2> ${base}_kraken.log
        gzip ${base}_kraken.txt
      """    
}

process braken2 {
    container "quay.io/biocontainers/bracken:2.6.1--py39h7cff6ad_2"
    containerOptions = "-v /bioinfonas/databases:/bioinfonas/databases:ro -v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
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
        bracken -d ${KRAKEN2_DB} -i ${base}_kraken.tsv -o ${base}_bracken_species.tsv -w ${base}_bracken_species_report.tsv -l S -t 10 > ${base}_bracken_species.log
        bracken -d ${KRAKEN2_DB} -i ${base}_kraken.tsv -o ${base}_bracken_genus.tsv -w ${base}_bracken_genus_report.tsv -l G -t 10 > ${base}_bracken_genus.log
        /scripts/create_import.py -ds ${md.ds} -o ${base}_import_taxa.csv -g ${base}_bracken_genus.tsv -sp ${base}_bracken_species.tsv -k ${kraken_report}
      """    
}

workflow step_3TX_class__kraken2 {
    take: reads
    main:
      kraken2(reads).report | braken2
     emit:
       genus_report = braken2.out.genus_report
}

workflow {
  step_3TX_class__kraken2(getSingleInput())
}
