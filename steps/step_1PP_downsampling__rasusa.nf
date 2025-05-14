nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata } from '../functions/common.nf'
include { getSingleInput;optionalOption;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '1PP_downsampling'
def METHOD = 'rasusa' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"
def PARAMS_PREFIX = "${STEP}__${METHOD}___"

process rasusa {
    container "quay.io/biocontainers/rasusa:2.0.0--h031d066_0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"   
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [opts:opts])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fastq.gz'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*_input.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_rasusa"

      opts = optionalOption('--num', METHOD + '__num')
      opts += optionalOption('--bases', METHOD + '__bases')
      opts += optionalOption('--genome-size', METHOD + '__genome_size')
      opts += optionalOption('--coverage', METHOD + '__coverage')
      opts += optionalOption('--frac', METHOD + '__frac') 

       if (r2) { 
        """
        rasusa reads ${opts} -s 1 ${r1} ${r2} -o ${base}_R1.fastq.gz -o ${base}_R2.fastq.gz -v
        """
      } else {
        """
        rasusa reads ${opts} -s 1 ${r1} -o ${base}.fastq.gz -v 
        """        
      }     
}

workflow step_1PP_downsampling__rasusa {
    take: 
      reads    
    main:
      rasusa(reads)
}

workflow {
    step_1PP_downsampling__rasusa(getSingleInput())
}