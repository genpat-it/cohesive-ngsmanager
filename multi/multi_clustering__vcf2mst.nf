nextflow.enable.dsl=2

include { taskMemory;filename } from '../functions/common.nf'
include { getVCFs;param;optional;optionalBoolean } from '../functions/parameters.nf'

def METADATA_ID_FULL = param('metadata_with_full_id') as boolean

process bcftools {
    container "quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      tuple val(riscd), path(vcf)
    output:
      tuple val(riscd), path("${vcf}.filtered"), emit: filtered_vcf
    script:
      """
      bcftools norm -a ${vcf} -o decomposed_complex.vcf
      bcftools view -v snps decomposed_complex.vcf -o snps_only.vcf
      bcftools norm -m +any snps_only.vcf -o "${vcf}.filtered"
	    """
}

process vcf2mst {
    container "${LOCAL_REGISTRY}/bioinfo/vcf2mst:0.0.1--d587d682e9"
    memory { taskMemory( 4.GB, task.attempt ) }
    input:
      path(summary)
      path(vcf_files)
    output:
      path '*'
      path 'HDmatrix.tsv', emit: matrix
      path '*.sh', hidden: true
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'tree.nwk'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "vcf2mst.cfg" }
    script:
      extra = optional('multi_clustering__vcf2mst__extra')
      """
        while IFS=\$'\t' read -r key file; do
            if (${METADATA_ID_FULL}); then
                KEY=\$key
            else
                KEY=\$(sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_.+/\\1/' <<< "\$file")
            fi
            mv \$file \${KEY}
            echo "\${KEY}" 
        done < ${summary} > vcf_files.txt

        vcf2mst.pl vcf_files.txt tree.nwk vcf ${extra} > vcf2mst.log
        cp /tmp/hamming_distance_matrix.tsv HDmatrix.tsv
      """
}

process dists {
    container "quay.io/biocontainers/cgmlst-dists:0.4.0--hec16e2b_2"
    memory { taskMemory( 5.GB, task.attempt ) }
    input:
      path(matrix)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{*.csv}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "cgmlst-dists.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "cgmlst-dists.cfg" }
    script:
      """
        cgmlst-dists -c ${matrix} > vcf2mst_dists_matrix.csv
      """
}


workflow multi_clustering__vcf2mst {
    take: 
      input
    main:
      if (optionalBoolean('multi_clustering__vcf2mst__vcf_normalization')) {
        vcf_input = bcftools(input).filtered_vcf
      } else {
        vcf_input = input
      }
      summary = vcf_input.collectFile { [ "summary.tsv", it[0] + '\t' + filename(it[1]) + '\n' ] }
      files= vcf_input.collect { it[1] }

      matrix=vcf2mst(summary, files).matrix
      dists(matrix)
}

workflow {
    multi_clustering__vcf2mst(getVCFs())
}
