nextflow.enable.dsl=2

include {  parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { param;getInput;optional } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '3TX_species'
def METHOD = 'blast' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def BLAST_DB_PATH = param('step_3TX_species__blast___blastdb')

def BLAST_OUTFMT = [
  "6" : [ ext : 'tsv', cols: 'qaccver saccver pident stitle sskingdom ssciname staxid length mismatch gapopen qstart qend sstart send evalue bitscore' ],
  "7" : [ ext : 'tsv', cols: 'qaccver saccver pident stitle sskingdom ssciname staxid length mismatch gapopen qstart qend sstart send evalue bitscore' ],
  "8" : [ ext : 'asn1.txt'],
  "9" : [ ext : 'asn1'],
  "10" : [ ext : 'csv', cols: 'qaccver saccver pident sskingdom ssciname staxid length mismatch gapopen qstart qend sstart send evalue bitscore' ],
  "11" : [ ext : 'blast.asn1'],
  "12" : [ ext : 'json'],
  "13" : [ ext : 'json'],
  "15" : [ ext : 'json'],
  "17" : [ ext : 'sam', cols: 'qaccver saccver pident stitle sskingdom ssciname staxid length mismatch gapopen qstart qend sstart send evalue bitscore' ]
]

process resolve_taxids {
    container 'ncbi/blast:2.14.0'
    containerOptions = " -v ${BLAST_DB_PATH}:/blast/blastdb:ro -u 0:0"
    memory { taskMemory(4.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(scaffolds200)
      val(taxid)
    output:
      path '*'
      tuple  val(taxid), path("${taxid}.txids"), emit: taxdata
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.txids'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_resolve_taxids"
      """    
        /blast/bin/get_species_taxids.sh -t ${taxid} > ${taxid}.txids  2> ${base}.log   
        [ -s ${taxid}.txids ]
      """
}

process blastn {
    container 'ncbi/blast:2.14.0'
    containerOptions = " -v ${BLAST_DB_PATH}:/blast/blastdb:ro -u 0:0 --network none "
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    cpus { [24, params.max_cpus as int].min() }
    stageInMode 'copy'
    when:
      db == 'nt'
    input:
      tuple val(riscd_input), path(scaffolds200)
      val(db)
      tuple val(taxid), path(taxids_file)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [db: db, taxid: taxid])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      max_target_seqs = param('step_3TX_species__blast___max_target_seqs')
      format = param('step_3TX_species__blast___outfmt')
      extra = optional('step_3TX_species__blast___extra')
      ext = BLAST_OUTFMT["${format}"]?.ext ?: 'result'
      cols = BLAST_OUTFMT["${format}"]?.cols ?: ''
      """    
        blastn -query ${scaffolds200} \
          -db ${db} \
          -num_threads ${task.cpus} \
          -taxidlist ${taxids_file} \
          -mt_mode 1 \
          -outfmt '${format} ${cols}' \
          -max_target_seqs ${max_target_seqs} \
          -out ${base}_${db}.${ext} \
          ${extra} &> ${base}.log
      """
}

workflow step_3TX_species__blast {
    take: 
      data
      db
      taxid
    main:
      taxdata = resolve_taxids(data, taxid).taxdata
      blastn(data, db, taxdata)
}

workflow {
    step_3TX_species__blast(getInput(), param('step_3TX_species__blast___db'), param('step_3TX_species__blast___taxid'))
}