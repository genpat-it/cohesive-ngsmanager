nextflow.enable.dsl=2

include { taskMemory;flattenPath } from '../functions/common.nf'
include { getReferenceUnkeyed;getResult;getInput;param } from '../functions/parameters.nf'

def GEO_RESOLUTION_COLUMNS = param('multi_clustering__reportree__summary_geo_column')
def SUMMARY_DATE_ALIASES = param('multi_clustering__reportree__summary_date_aliases')
def SAMPLE_COLUMN = param('multi_clustering__reportree__summary_sample_column')

IMAGES = [
  '2.2.1': 'staphb/cfsan-snp-pipeline:2.2.1',
  '2.0.2': 'cfsanbiostatistics/snp-pipeline@sha256:448787923371ade95217982814db25efb1e01287a8180b523d76a9f093f97d01'
]

def DOCKER_IMAGE = IMAGES[param('multi_clustering__cfsan__version')] ?: (exit 2, "params (multi_clustering__cfsan__version) not valid");

process cfsan_snp_pipeline {
    container DOCKER_IMAGE
    containerOptions = "-v ${workflow.projectDir}/scripts/multi_clustering__cfsan:/scripts:ro"
    input:
      path(samples)
      tuple val(_), val(reference), path(refPath)
    output:
      path '**'
      path 'results/snpma.fasta', emit: snpma
      path '{*.sh,*.log}', hidden: true
    stageInMode 'symlink'  
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'results/*.tsv', saveAs: { filename -> flattenPath(filename) }      
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'results/*.vcf', saveAs: { filename -> flattenPath(filename) }      
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'results/*.fasta', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "cfsan_snp_pipeline.cfg" }
    script:
      version = param('multi_clustering__cfsan__version')
      """
        trap "find  \\( -name '*.bam' -o -name '*.sam' -o -name '*.pileup' \\) -delete" EXIT
        for FILE in ${samples}; do CMP=`echo \$FILE | sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_[[:graph:]]+/\\1/'`; mkdir -p samples/\${CMP} && mv \$FILE samples/\${CMP}/ ; done
        cfsan_snp_pipeline run -c /scripts/snppipeline_${version}.conf -m soft -o results --samples_dir samples ${refPath} >> cfsan_snp_pipeline.log
      """
}

process iqtree {
    container "quay.io/biocontainers/iqtree:1.6.12--he513fc3_0"
    stageInMode 'copy'
    input:
      path(snpma)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path("snpma.fasta.treefile"), emit: nwk
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'snpma.fasta.treefile', saveAs: { "snpma.fasta.nwk" } 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "iqtree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "iqtree.cfg" }
    script:
      """
        iqtree -s ${snpma} -nt AUTO
      """
}

process augur {
    container "quay.io/biocontainers/augur:22.0.0--pyhdfd78af_0"
    memory { taskMemory( 4.GB, task.attempt ) }
    input:
      path(nwk)
      path(metadata)
      path(geodata)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'auspice.json'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "augur.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "augur.cfg" }
    script:
      """
        cat ${metadata} | sed 's/${SAMPLE_COLUMN}/name/i' \
           | sed 's/${SUMMARY_DATE_ALIASES}/date/i' > augur_metadata.tsv
        METADATA_LIST=\$(head -n 1 augur_metadata.tsv | tr \$'\t' ' ')
        augur refine --tree ${nwk} --output-tree tree_tt.nwk --output-node-data refine.node.json --metadata augur_metadata.tsv
        augur export v2 --tree tree_tt.nwk --node-data refine.node.json --output auspice.json \
          --color-by-metadata \${METADATA_LIST} \
          --geo-resolutions ${GEO_RESOLUTION_COLUMNS} \
          --metadata augur_metadata.tsv \
          --lat-longs ${geodata}
      """
}

workflow multi_clustering__cfsan {
    take: 
        input
        reference
        metadata
        geodata
    main:
        snpma = cfsan_snp_pipeline(input, reference).snpma
        nwk = iqtree(snpma).nwk
        augur(nwk, metadata, geodata)
}

workflow {
    reads = getInput()
        .map { it[1] }
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatten()     
        .collect()  
    multi_clustering__cfsan(reads, getReferenceUnkeyed('fa'), param('metadata'), param('geodata'));
}