nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata } from '../functions/common.nf'
include { getSingleInput;param;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '1PP_downsampling'
def METHOD = 'bbnorm' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process bbnorm {
    container "quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      isCompatibleWithSeqType(reads, 'illumina_paired', task.process)        
    input:
      tuple val(riscd_input), path(reads)
      val(k)
      val(target)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [k:k, target:target])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fastq.gz'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.hist,*_input.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = reads
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_bbnorm_k${k}_t${target}"
      // bbnorm has trouble calculating the max heap size inside a container
      javaHeapLimit = ((params.max_memory as nextflow.util.MemoryUnit).getMega() * 0.85) as int
      """
        bbnorm.sh \
          in=${r1} \
          in2=${r2} \
          out=${base}_R1.fastq.gz \
          out2=${base}_R2.fastq.gz \
          hist=${base}.hist \
          k=${k} \
          target=${target} \
          -Xmx${javaHeapLimit}m
      """
}

workflow step_1PP_downsampling__bbnorm {
    take: 
      reads
      k
      target
    main:
      bbnorm(reads, k, target)
}

workflow {
    step_1PP_downsampling__bbnorm(getSingleInput(), param('k'), param('target'))
}