nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;taskTime } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '2AS_denovo'
def METHOD = 'shovill' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process shovill {
    container "${LOCAL_REGISTRY}/bioinfo/shovill:1.1.0--d84470570e"    
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 16.GB, task.attempt ) }
    time { taskTime( 45.m, task.attempt ) }    
    cpus 8   
    when:
      isCompatibleWithSeqType(reads, ['illumina_paired'], null)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '**'
      tuple val(riscd), path("${base}.fasta"), emit: assembly
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    script:
      (t1,t2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
        shovill --outdir out --minlen 200 --cpus ${task.cpus} --ram ${task.memory.toGiga()} --R1 ${t1} --R2 ${t2}
        mv out/contigs.fa ${base}.fasta
      """
}

process shovill_se {
    container "${LOCAL_REGISTRY}/bioinfo/shovill-se:1.1.1--ba51ea69e5"    
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 16.GB, task.attempt ) }
    time { taskTime( 45.m, task.attempt ) }    
    cpus 8   
    when:
      isCompatibleWithSeqType(reads, ['ion'], null)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '**'
      tuple val(riscd), path("${base}.fasta"), emit: assembly
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    script:
      (t1,t2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
        shovill-se --outdir out --minlen 200 --cpus ${task.cpus} --ram ${task.memory.toGiga()} --se ${t1} --opts '--sc --iontorrent' --kmers '31,33,55'
        mv out/contigs.fa ${base}.fasta
      """        
}

process quast {
    container 'quay.io/biocontainers/quast:4.4--boost1.61_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    time { taskTime( 5.m, task.attempt ) }    
    input:
      tuple val(_), path(assembly)
    output:
      path '*_quast.*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_quast"
      """
      quast -m 200 --fast -o quast ${assembly} > ${base}.log ;
			cut -f1,14,15,16,17,18,19,20,21 quast/transposed_report.tsv > ${base}.csv ; 
      """
}

process checkm {
    container 'quay.io/biocontainers/checkm-genome:1.1.3--py_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 48.GB, task.attempt ) }
    time { taskTime( 15.m, task.attempt ) }    
    cpus 16
    input:
      tuple val(_), path(assembly)
    output:
      path '{out/**,*.tab}'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: 'results.tab', saveAs: { "${base}_results.tab" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: 'out/files/storage/bin_stats.analyze.tsv', saveAs: { "${base}_bin_stats.analyze.tsv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: 'out/files/checkm.log', saveAs: { "${base}_hmm.log" }         
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_checkm"
      """
        mkdir in && mkdir -p out/files && cp ${assembly} in/
        checkm lineage_wf -x fasta -t ${task.cpus} --reduced_tree --tab_table ./in ./out/files | tee results.tab.tmp
        grep "^Bin" -A1 results.tab.tmp > results.tab || true && rm results.tab.tmp
      """
}

workflow step_2AS_denovo__shovill {
    take: rawreads
    main:
      contigs_from_pe = shovill(rawreads).assembly
      contigs_from_se = shovill_se(rawreads).assembly
      contigs = contigs_from_pe.mix(contigs_from_se)
      quast(contigs)          
      if (!params.skip_checkm) {
        checkm(contigs)
      }
    emit:
      assembly = contigs
}

workflow {
    step_2AS_denovo__shovill(getSingleInput())
}
