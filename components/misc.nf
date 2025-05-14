include { taskMemory } from '../functions/common.nf'

process sequence_merge {
    container "ubuntu:20.04"
    memory { taskMemory( 1.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'sequence_merge.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "sequence_merge.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "sequence_merge.cfg" }
    input:
      path(sequence1)
      path(sequence2)
    output:
      path '*'
      path 'sequence_merge.fasta', emit: merged_sequences
      path '{*.sh,*.log}', hidden: true
    script:
      """
        awk 1 ${sequence1} ${sequence2} > sequence_merge.fasta
      """
}

process metadata_merge {
    container "quay.io/biocontainers/csvtk:0.31.0--h9ee0642_0"
    memory { taskMemory( 1.GB, task.attempt ) }
    cpus 4
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'metadata_all.tsv'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "metadata_merge.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "metadata_merge.cfg" }
    input:
      path(metadata)
      path(metadata_ext)
    output:
      path '*'
      path 'metadata_all.tsv', emit: merged_metadata
      path '{*.sh,*.log}', hidden: true
    script:
      """
        csvtk concat -t -i ${metadata} ${metadata_ext} > metadata_all.tsv 
      """
}

process geodata_merge {
    container "ubuntu:20.04"
    memory { taskMemory( 1.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'geodata_merge.tsv'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "geodata_merge.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "geodata_merge.cfg" }
    input:
      path(f1)
      path(f2)
    output:
      path '*'
      path 'geodata_merge.tsv', emit: merged_geodata
      path '{*.sh,*.log}', hidden: true
    script:
      """
        awk 1 ${f1} ${f2} > geodata_merge.tsv
      """
}