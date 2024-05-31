nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '2AS_denovo'
def METHOD = 'unicycler' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process unicycler {
    container 'docker.io/biocontainers/unicycler:v0.4.7dfsg-2-deb_cv1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 20.GB, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, 'illumina_paired', task.process)    
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path("${base}_unicycler_assembly.fasta"), emit: scaffolds
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "{*.fasta,*.gfa}"    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_spades.cfg" }
    script:
      (t1,t2) = reads
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      unicycler -1 ${t1} -2 ${t2} -o . -t 16 > ${base}_unicycler.log 2>&1
      mv assembly.fasta ${base}_unicycler_assembly.fasta
      """
}

process assembly_filter {
    container "ghcr.io/genpat-it/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    memory { taskMemory( 3.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      path(scaffolds)
    output:
      tuple val(riscd), path("${base}_${METHOD}_scaffolds_L200.fasta"), emit: fasta
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_assemblyfilter.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
      /scripts/AssemblyFilter.py -n ${base} -f ${scaffolds} -u -l 200 -c 0 ;
      """
}

process quast {
    container 'quay.io/biocontainers/quast:4.4--boost1.61_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      tuple val(_), path(l200)
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

workflow step_2AS_denovo__unicycler {
    take: data
    main:
      unicycler(data).scaffolds
      assembly_filter(unicycler.out.scaffolds).fasta | quast
    emit:
      assembled = assembly_filter.out.fasta       
}

workflow {
    step_2AS_denovo__unicycler(getSingleInput())
}

