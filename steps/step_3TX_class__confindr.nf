nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;isSpeciesSupported;taskMemory } from '../functions/common.nf'
include { getSingleInput;param;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '3TX_class'
def METHOD = 'confindr' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_ALLOWED = [
  'escherichia',
  'salmonella',
  'listeria' 
]

process confindr {
    container "${LOCAL_REGISTRY}/bioinfo/confindr:0.7.4--8fe19dd246"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 24.GB, task.attempt ) }
    cpus  16
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process) && isSpeciesSupported(genus_species, GENUS_ALLOWED, reads, task.process)
    input:
      tuple val(riscd_input), path(reads)
      val(genus_species)
    output:
      path '*'
      path '*_report.csv'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'confindr_report.csv', saveAs: { "${base}_report.csv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*_rmlst.csv,*_contamination.csv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: 'confindr_log.txt', saveAs: { "${base}_log.txt" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_confindr"     
      """
        confindr.py \
          -i ./ \
          -o . \
          -d /db \
          -t ${task.cpus} \
          --cross_details
      """   
}

workflow step_3TX_class__confindr {
    take: 
      reads
      genus_species
    main:
      confindr(reads, genus_species)
}

workflow {
    step_3TX_class__confindr(getSingleInput(),param('genus_species'))
}