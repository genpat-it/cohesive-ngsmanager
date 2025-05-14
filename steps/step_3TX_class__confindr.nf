nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;isSpeciesSupported;taskMemory } from '../functions/common.nf'
include { getSingleInput;param;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '3TX_class'
def METHOD = 'confindr' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"


process confindr {
    container "${LOCAL_REGISTRY}/bioinfo/confindr:0.7.4--9508ee4865"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 24.GB, task.attempt ) }
    cpus  16
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'confindr_report.csv', saveAs: { "${base}_report.csv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*_rmlst.csv', saveAs: { "${base}_rmlst.csv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*_contamination.csv', saveAs: { "${base}_contamination.csv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: 'confindr_log.txt', saveAs: { "${base}_log.txt" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"     
      if (r2) { 
        """
        mkdir input && cd input
        ln -s ../$r1 ${md.cmp}_R1.fastq.gz 
        ln -s ../$r2 ${md.cmp}_R2.fastq.gz       
        cd ..
        confindr.py \
          -i ./input \
          -o . \
          -d /db \
          -t ${task.cpus} \
          --cross_details
        """
      } else {
        """
        mv $r1 ${md.cmp}_R1.fastq.gz 
        confindr.py \
          -i ./ \
          -o . \
          -d /db \
          -t ${task.cpus} \
          --cross_details
        """        
      }                
}

workflow step_3TX_class__confindr {
    take: 
      reads
    main:
      confindr(reads)
}

workflow {
    step_3TX_class__confindr(getSingleInput())
}