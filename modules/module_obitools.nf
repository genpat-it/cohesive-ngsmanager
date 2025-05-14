nextflow.enable.dsl=2

include { taskMemory;parseMetadataFromFileName;getRisCd;stepInputs;executionMetadata;extractKey } from '../functions/common.nf'
include { getInput;param;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()

def DB_PATH = param('module_obitools__db_path')
def DB_RAW = "${DB_PATH}/${param('module_obitools__db')}"

process trimmomatic {
    container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      isCompatibleWithSeqType(reads, ['illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      tuple val(riscd_input), path("${base}_R*_trimmomatic.fastq.gz"), emit: fastq
      path '{*trimmomatic.log,*.json}'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fastq.gz'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}_trimmomatic.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}_trimmomatic.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      min_len = param('module_obitools__trimmomatic_min_len')  
      sliding_window = param('module_obitools__trimmomatic_sliding_window')  
      """
        trimmomatic PE -threads 2 $r1 $r2 ${base}_R1_trimmomatic.fastq.gz ${base}_R1_unpaired.fastq.gz ${base}_R2_trimmomatic.fastq.gz ${base}_R2_unpaired.fastq.gz \
        SLIDINGWINDOW:${sliding_window} \
        MINLEN:${min_len} 2>> ${base}_trimmomatic.log;
      """
          
}
/*
        ILLUMINACLIP:/usr/local/share/trimmomatic-0.36-6/adapters/NexteraPE-PE.fa:2:30:10 \
                LEADING:optional('module_obitools__')   \
                        TRAILING:25 \

*/

//FASTQC forse dopo
process illuminapairedend {
    container "quay.io/biocontainers/obitools:1.2.13--py27h9801fc8_4"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}_illuminapairedend.fastq"), emit: fastq
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fastq'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}_illuminapairedend.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}_illuminapairedend.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      score_min = param('module_obitools__illuminapairedend_score_min')  
      """
        illuminapairedend --sanger --without-progress-bar --score-min='${score_min}' -r ${r2} ${r1} > ${base}_illuminapairedend.fastq
      """
}

process obigrep {
    container "quay.io/biocontainers/obitools:1.2.13--py27h9801fc8_4"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}_obigrep.fastq"), emit: fastq
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fastq'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}_obigrep.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}_obigrep.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
         obigrep --sanger -p 'mode!="joined"' ${reads}  > ${base}_obigrep.fastq
      """
}

process obiuniq {
    container "quay.io/biocontainers/obitools:1.2.13--py27h9801fc8_4"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}_obiuniq.fasta"), emit: fasta
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}_obiuniq.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}_obiuniq.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
        obiuniq --without-progress-bar -m 'sample' --sanger  ${reads} > ${base}_obiuniq.fasta
      """
}

process obiannotate {
    container "quay.io/biocontainers/obitools:1.2.13--py27h9801fc8_4"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}_obiannotate.fasta"), emit: fasta
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}_obiannotate.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}_obiannotate.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
        obiannotate --without-progress-bar --keep='count' --keep='merged_sample' --length  --fasta ${reads} > ${base}_obiannotate.fasta
      """
}


process obigrep2 {
    container "quay.io/biocontainers/obitools:1.2.13--py27h9801fc8_4"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}_obigrep2.fasta"), emit: fasta
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}_obigrep2.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}_obigrep2.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
         obigrep --without-progress-bar  -p count\\>=10 --fasta ${reads} > ${base}_obigrep2.fasta
      """
}

process obiclean {
    container "quay.io/biocontainers/obitools:1.2.13--py27h9801fc8_4"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}_obiclean.fasta"), emit: fasta
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}_obiclean.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}_obiclean.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      d = param('module_obitools__obiclean_d')  
      r = param('module_obitools__obiclean_r')  
      """
        obiclean --without-progress-bar -d '${d}' -r '${r}' -H --fasta ${reads} > ${base}_obiclean.fasta
      """
}

process blast {
    container "ncbi/blast:2.15.0"
    containerOptions = " -v ${params.module_obitools__db_path}:/blast/blastdb:ro -u 0:0 --network none "
    memory { taskMemory( 5.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    cpus 16
    input:
      tuple val(riscd_input), path(reads)
      val(db_name)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}_blast.tsv"), emit: tsv
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.tsv'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}_blast.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}_blast.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      evalue = param('module_obitools__blast_evalue')  
      identity = param('module_obitools__blast_identity')  
      db_name = param('module_obitools__db')  
      """
        blastn -query ${reads} -db ${db_name} -task 'megablast' -evalue '${evalue}' -out  ${base}_blast.tsv -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles'  -num_threads "${task.cpus}" -strand both -dust yes  -max_target_seqs '1' -perc_identity '${identity}'
      """
}

process db_match_id {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/module_obitools:/scripts:ro"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(blast_result)
      path(db_raw)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path '*.sh', hidden: true
      tuple val(riscd_input), path("${base}.fasta"), emit: fasta
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(blast_result.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_db_match"
      """
        /scripts/seq_filter_by_id.py -i ${db_raw} -f 'fasta' -p ${base}.fasta  -l UNION ${blast_result} '2'
	    """
}


process db_match_id2 {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/module_obitools:/scripts:ro"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(blast_result)
      tuple val(_), path(fasta_clean)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}.fasta"), emit: fasta     
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(blast_result.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_fasta_match"
      """
        /scripts/seq_filter_by_id.py -i ${fasta_clean} -f 'fasta' -p ${base}.fasta  -l UNION ${blast_result} '1'
	    """
}

process obitab {
    container "quay.io/biocontainers/obitools:1.2.13--py27h9801fc8_4"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(fasta)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}.tsv"), emit: tsv
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.tsv'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(fasta.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_db_match"
      """
        obitab --without-progress-bar --fasta  ${fasta} > ${base}.tsv
      """
}

process obitab2 {
    container "quay.io/biocontainers/obitools:1.2.13--py27h9801fc8_4"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(fasta)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}.tsv"), emit: tsv
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.tsv'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(fasta.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_fasta_match"
      """
        obitab --without-progress-bar --fasta  ${fasta} > ${base}.tsv
      """
}

process join_tabs {
    container "bgruening/galaxy-stable"
    containerOptions = "-u 0:0"
    memory { taskMemory( 2.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(tsv1)
      tuple val(_), path(tsv2)
      tuple val(__), path(blast_result)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      tuple val(riscd_input), path("${base}.tsv"), emit: fasta
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.tsv'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(tsv1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_join_tabs"
      """
        #!/bin/bash -euo pipefail -l
        conda activate base 
        python /galaxy-central/tools/filters/join.py ${tsv2} ${blast_result} 1 1 tmp.file \
        --index_depth=3 --buffer=50000000 
        python /galaxy-central/tools/filters/join.py tmp.file ${tsv1} 25 1 ${base}.tsv \
        --index_depth=3 --buffer=50000000 
      """
}

workflow module_obitools {
    take: 
      raw_reads
      blastdb
      db_raw
    main:  
      trimmed = trimmomatic(raw_reads).fastq
      fastq_paired = illuminapairedend(trimmed).fastq
      fastq_grep = obigrep(fastq_paired).fastq
      fasta_uniq = obiuniq(fastq_grep).fasta
      fasta_annotated = obiannotate(fasta_uniq).fasta
      fasta_grep2 = obigrep2(fasta_annotated).fasta
      fasta_clean = obiclean(fasta_grep2).fasta
      blast_result = blast(fasta_clean, blastdb).tsv
      db_match = db_match_id(blast_result, db_raw).fasta
      fasta_match = db_match_id2(blast_result, fasta_clean).fasta
      tsv_db_batch = obitab(db_match).tsv
      tsv_fasta_batch = obitab2(fasta_match).tsv
      join_tabs(tsv_db_batch,tsv_fasta_batch,blast_result)     
  }

workflow {
    raw_reads = getInput()
    module_obitools(raw_reads, DB_PATH, DB_RAW)
}
