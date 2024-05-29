nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;extractKey } from '../functions/common'
include { getTrimmedReads;getParamTaxaId;getParamIncludeChildren;getParamIncludeParents;getKrakenResults } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '1PP_filtering'
def METHOD = 'krakentools' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process krakentools {
    container 'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0'
    memory { taskMemory( 1.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(krakenResults)
      tuple val(riscd_input2), path(reads)
      val(taxaid)
      val(include_parents)
      val(include_children)
    output:
      path '*'
      path '{*.sh,*.log,*.err}', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_input2], md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fastq.gz'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.err', saveAs: { "${base}_krakentools.err" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}_krakentools.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_krakentools.cfg" }
    script:
      (t1,t2) = reads
      (krakenReport,krakenOutput) = krakenResults
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      options = (include_parents ? ' --include-parents ' : '') + (include_children ? ' --include-children ' : '')
      suffix = "t${taxaid}" + (include_parents ? 'p' : '') + (include_children ? 'c' : '')
      """
        zcat ${krakenOutput} > kraken_output
        extract_kraken_reads.py --fastq-output -k kraken_output ${options} -r ${krakenReport} -1 ${t1} -2 ${t2} -t ${taxaid} -o ${base}_R1_krakentools_${suffix}.fastq -o2 ${base}_R2_krakentools_${suffix}.fastq
        gzip *${suffix}.fastq && rm kraken_output
      """
}


workflow step_1PP_filtering__krakentools {
    take: 
      kraken
      trimmed
      taxaid
      include_children
      include_parents
    main:
      krakentools(kraken, trimmed, taxaid, include_children, include_parents)
}

workflow {
    getKrakenResults().cross(getTrimmedReads(false)) { extractKey(it) }.multiMap { 
      kraken: it[0]
      trimmed: it[1]
    }.set { krakenAndTrimmed }
    step_1PP_filtering__krakentools(krakenAndTrimmed.kraken, krakenAndTrimmed.trimmed, getParamTaxaId(), getParamIncludeChildren(), getParamIncludeParents())
}