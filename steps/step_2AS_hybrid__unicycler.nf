nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;getLongReads;param } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '2AS_hybrid'
def METHOD = 'unicycler' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process unicycler {
    container 'quay.io/biocontainers/unicycler:0.5.0--py38h5cf8b27_3'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    cpus params.max_cpus
    when:
      isCompatibleWithSeqType(short_reads, 'illumina_paired', "${ENTRYPOINT} short_reads") // not robust enough && isCompatibleWithSeqType(long_reads, 'nanopore', "${ENTRYPOINT} long_reads")  
    input:
      tuple val(riscd_input), path(short_reads)
      tuple val(riscd_input2), path(long_reads)
    output:
      path '*'
      path("${base}_assembly.fasta"), emit: scaffolds
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_input2], md, ex, STEP, METHOD, [long_reads:riscd_input2])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "*.fasta"    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "assembly.gfa", saveAs: { "${base}.gfa" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (t1,t2) = short_reads
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      mode = param('step_2AS_hybrid__unicycler__mode')      
      """
      unicycler -1 ${t1} -2 ${t2} -l ${long_reads} -o . -t ${task.cpus} --mode ${mode}
      mv assembly.fasta ${base}_assembly.fasta
      """
}


process quast {
    container 'quay.io/biocontainers/quast:4.4--boost1.61_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      path(l200)
    output:
      path '*_quast.*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_quast.cfg" }
    script:
      md = parseMetadataFromFileName(l200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      quast -m 200 --fast -o quast ${l200} > ${base}_quast.log ;
			cut -f1,14,15,16,17,18,19,20,21 quast/transposed_report.tsv > ${base}_quast.csv ; 
      """
}

workflow step_2AS_hybrid__unicycler {
    take: 
      short_reads
      long_reads
    main:
      scaffolds = unicycler(short_reads, long_reads).scaffolds
      quast(scaffolds)
    emit:
      scaffolds = scaffolds
}

workflow {  
    step_2AS_hybrid__unicycler(getSingleInput(), getLongReads())
}

