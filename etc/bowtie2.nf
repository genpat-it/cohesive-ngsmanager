nextflow.enable.dsl=2

include { logHeader } from '../functions/common.nf'
include { param } from '../functions/parameters.nf'

log.info logHeader('NGSMANAGER')

process bowtie2_index {
    container "ghcr.io/genpat-it/bowtie2:2.1.0--37ad014737"
    input:
      val(ref)
      path(fasta)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${ref}/fasta/result", pattern: '*.bt*'    
    publishDir mode: 'rellink', "${params.outdir}/${ref}/fasta/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${ref}/fasta/meta", pattern: '.command.sh', saveAs: { "${ref}.cfg" }
    script:
      base = "bowtie2_${ref}"
      """
      bowtie2-build ${fasta} ${ref} &> ${base}.log
      """
}

workflow build_index {
  bowtie2_index(param('reference'), param('reference_path'))
}