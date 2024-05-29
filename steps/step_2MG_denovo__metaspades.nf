nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '2MG_denovo'
def METHOD = 'metaspades' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process metaspades {
    container 'quay.io/biocontainers/spades:3.11.1--py27_zlib1.2.8_0'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 24.GB, task.attempt ) }
    cpus { [16, params.max_cpus as int].min() }       
    when:
      isCompatibleWithSeqType(reads, 'illumina_paired', task.process)    
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path("${base}_spades_scaffolds.fasta"), emit: scaffolds
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_spades.cfg" }
    script:
      (t1,t2) = reads
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      spades.py --meta -k 21,33,55,77 -t ${task.cpus}  -o spades -1 ${t1} -2 ${t2} > ${base}_spades.log 2>&1 
			mv spades/scaffolds.fasta ${base}_spades_scaffolds.fasta
      mv spades/contigs.fasta ${base}_spades_contigs.fasta
      """
}

process assembly_filter {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    memory { taskMemory( 3.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      path(scaffolds)
    output:
      tuple val(riscd), path("${base}_spades_scaffolds_L200.fasta"), emit: fasta
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_assemblyfilter.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
      /scripts/AssemblyFilter.py -n ${base} -f ${scaffolds} -l 200 -c 0 ;
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

workflow step_2MG_denovo__metaspades{
    take: data
    main:
      metaspades(data)
      assembly_filter(metaspades.out.scaffolds).fasta | quast
    emit:
      assembled = assembly_filter.out.fasta
}

workflow {
    step_2MG_denovo__metaspades(getSingleInput())
}

