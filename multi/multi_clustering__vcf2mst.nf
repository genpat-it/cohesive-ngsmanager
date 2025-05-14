nextflow.enable.dsl=2

include { taskMemory } from '../functions/common.nf'
include { getVCFs } from '../functions/parameters.nf'

process vcf2mst {
    container "ghcr.io/genpat-it/vcf2mst:0.0.1--d587d682e9"
    memory { taskMemory( 4.GB, task.attempt ) }
    input:
      path(vcf_files)
    output:
      path '*'
      path 'HDmatrix.tsv', emit: matrix
      path '*.sh', hidden: true
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'tree.nwk'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "vcf2mst.cfg" }
    script:
      """
        mkdir vcf_files 
        for FILE in ${vcf_files}; do CMP=`echo \$FILE | sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([[:digit:]]+\\.[[:alnum:]]+\\.[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+).+/\\1/'`; cp \$FILE vcf_files/\${CMP} ; done
        find vcf_files -mindepth 1 > vcf_list            
        vcf2mst.pl vcf_list tree.nwk vcf > vcf2mst.log
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
        matrix=vcf2mst(input).matrix
        dists(matrix)
}

workflow {
    getVCFs()
    .map { it[1] }
    .collect()
    .set { inputSet }
    multi_clustering__vcf2mst(inputSet)
}
