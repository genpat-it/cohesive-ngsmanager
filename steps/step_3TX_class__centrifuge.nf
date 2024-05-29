nextflow.enable.dsl=2

include { getEmpty;flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { param;isFullOutput;getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters'
include { stepInputs;getRisCd } from '../functions/common.nf'

def DB_PATH=param('step_3TX_class__centrifuge__db_path')
def DB_NAME=param('step_3TX_class__centrifuge__db_name')

def ex = executionMetadata()

def STEP = '3TX_class'
def METHOD = 'centrifuge' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process centrifuge {
    container "quay.io/biocontainers/centrifuge:1.0.4_beta--h9a82719_6"
    containerOptions = "-v ${DB_PATH}:${DB_PATH}:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    cpus { [16, params.max_cpus as int].min() }
    memory { taskMemory( 16.GB, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, ['nanopore'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      tuple path("${base}_genus.report"), path("${base}_species.report"), path("${base}_KrakenLike.report"), emit: reports
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.output,*.report}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(reads.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}" 
      """
        centrifuge -p ${task.cpus} -x ${DB_PATH}/${DB_NAME} -U ${reads} -S ${base}.output --report-file ${base}.report
        head -n 1 ${base}.report > ${base}_genus.report; grep -w "genus" ${base}.report | sort -t "\t" -nr -k5 >> ${base}_genus.report
        head -n 1 ${base}.report > ${base}_species.report; grep -w "species" ${base}.report | sort -t "\t" -nr -k5 >> ${base}_species.report         
        centrifuge-kreport -X ${DB_PATH}/${DB_NAME} ${base}.report > ${base}_KrakenLike.report
      """    
}

process import_taxa {
    container "quay.io/biocontainers/biopython:1.78"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    memory { taskMemory( 1.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple path(genus), path(species), path(krakenLike)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*_import_taxa.csv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_import_taxa.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}_import_taxa.log" }
    script:
      md = parseMetadataFromFileName(genus.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      /scripts/create_import.py -ds ${md.ds} -o ${base}_import_taxa.csv -g ${genus} -sp ${species} -k ${krakenLike}
      """
}

workflow step_3TX_class__centrifuge {
    take: reads
    main: 
      centrifuge(reads).reports | import_taxa
}

workflow {
  step_3TX_class__centrifuge(getSingleInput())
}
