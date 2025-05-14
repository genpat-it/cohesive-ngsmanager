nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent;optional } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '2AS_denovo'
def METHOD = 'flye' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process flye {
    container "quay.io/biocontainers/flye:2.9.2--py310h2b6aa90_2"    
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 32.GB, task.attempt ) }
    cpus 16   
    when:
      isCompatibleWithSeqType(reads, ['nanopore'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '**'
      tuple val(riscd), path("${base}.fasta"), emit: assembly
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: 'assembly_info.txt', saveAs: { "${base}_assembly_info.txt" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    script:
      (t1,t2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      g = "${optional('step_2AS_denovo__flye__genome_size')}".replaceAll(/[\s\n\t\r"'$\{]/, "")
      g_par = g ? " --genome-size ${g}" : ''      
      m_par = optional('step_2AS_denovo__flye__meta') ? " --meta " : ''          
      """
        flye \
        --nano-hq \
        ${reads} \
        --out-dir .\
        --threads ${task.cpus} ${g_par} ${m_par}
        mv assembly.fasta ${base}.fasta
      """        
}

process quast {
    container 'quay.io/biocontainers/quast:4.4--boost1.61_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      tuple val(_), path(assembly)
    output:
      path '*_quast.*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_quast"
      """
      quast -m 200 --fast -o quast ${assembly} ;
			cut -f1,14,15,16,17,18,19,20,21 quast/transposed_report.tsv > ${base}.csv ; 
      """
}


workflow step_2AS_denovo__flye {
    take: reads
    main:
      contigs = flye(reads).assembly
      quast(contigs)
    emit:
      assembly = contigs
}

workflow {
    step_2AS_denovo__flye(getSingleInput())
}