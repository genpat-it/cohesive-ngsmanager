nextflow.enable.dsl=2

include { getInput;param;optional } from '../functions/parameters.nf'
include { flattenPath } from '../functions/common.nf'

process panaroo {
    container "quay.io/biocontainers/panaroo:1.3.3--pyhdfd78af_0"
    input:
      path gffs
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'results/*', saveAs: { filename -> flattenPath(filename) }      
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "panaroo.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "panaroo.log" }
    script:
      """          
         panaroo -i *.gff \
            -o results \
            --clean-mode ${param('multi_pangenome__panaroo__clean_mode')} \
            --remove-invalid-genes \
            --threshold ${param('multi_pangenome__panaroo__threshold')} \
            --family_threshold ${param('multi_pangenome__panaroo__family_threshold')} \
            --len_dif_percent ${param('multi_pangenome__panaroo__len_dif_percent')} \
            -t ${param('multi_pangenome__panaroo__threads')} \
            --alignment core \
            --aligner mafft \
            --merge_paralogs  \
            ${optional('multi_pangenome__panaroo__extra')}
      """
}

workflow multi_pangenome__panaroo {
    take: 
        input
    main:      
        input
            .map { it[1] }
            .collect()
            .set { gffs }   
        panaroo(gffs)
}

workflow {
    multi_pangenome__panaroo(getInput())
}