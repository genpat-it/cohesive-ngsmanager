nextflow.enable.dsl=2

include { taskMemory;getEmpty } from '../functions/common.nf'
include { _getAlleles;param;optional;optionalOrDefault;getVCFs } from '../functions/parameters.nf'

def SUMMARY_DATE_ALIASES = param('multi_clustering__reportree__summary_date_aliases')  
def SUMMARY_COLUMNS = param('multi_clustering__reportree__summary_columns')  
def SAMPLE_COLUMN = param('multi_clustering__reportree__summary_sample_column')  
def GEO_RESOLUTION_COLUMNS = param('multi_clustering__reportree__summary_geo_column')  

process prepare_metadata {
    container "ubuntu:20.04"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      path(metadata)
    output:
      path '*'
      path 'reportree_metadata.tsv', emit: metadata
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'reportree_metadata.tsv' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "prepare_metadata.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "prepare_metadata.cfg" }
    script:
      """
         sed -E 's/${SUMMARY_DATE_ALIASES}/date/i' ${metadata} > reportree_metadata.tsv
      """
}

process maaft {
    container "quay.io/biocontainers/mafft:7.520--h031d066_3"
    cpus { [32, params.max_cpus as int].min() }   
    input:
      path(sequences)
    output:
      path '*'
      path 'alignment.fasta', emit: alignment
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{*.fasta}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "maaft.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "maaft.cfg" }
    script:     
      """
        for f in ${sequences}; do 
          cmp=`echo -n \$f | sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_[[:graph:]]+/\\1/'`
          awk "NR==1 {\\\$0=\\">\$cmp\\"}1" \$f >> input.fasta
        done   
        mafft \
        --thread ${task.cpus} \
        --auto \
        input.fasta > alignment.fasta
      """
}

process reportree_gt {
    container "ghcr.io/genpat-it/reportree:2.4.1--088b6651b8"
    cpus { [32, params.max_cpus as int].min() }   
    input:
      path(alignment)
      path(metadata)
      path(nomenclature_path), stageAs: 'prev_nomenclature.tsv'
    output:
      path '*'
      path 'gt_dist_hamming.tsv', emit: matrix
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'gt_zooms.txt', saveAs: { "zooms.txt" }
    publishDir mode: 'copy', "${params.outdir}/result", pattern: '{*.tsv,*.txt,*.nwk}'
    publishDir mode: 'copy', "${params.outdir}/result", pattern: '*_*[0-9]'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "reportree_gt.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "reportree_gt.cfg" }
    script:
      thr = param('multi_clustering__reportree__thr')      
      soi = optional('multi_clustering__reportree__sample_of_interest').replaceAll(/[\s\n\t\r"'$\{]/, "")
      soi_par = soi ? " --sample_of_interest ${soi}" : ''
      zoom = optional('multi_clustering__reportree__zoom_cluster_of_interest').replaceAll(/[\s\n\t\r"'$\{]/, "")
      zoom_par = zoom ? " --zoom-cluster-of-interest  ${zoom}" : ''
      zoom_subtree = optional('multi_clustering__reportree__subtree-of-interest').replaceAll(/[\s\n\t\r"'$\{]/, "")
      zoom_subtree_par = zoom_subtree ? " --subtree-of-interest  ${zoom_subtree}" : ''
      lociCalled = param('multi_clustering__reportree__loci_called')   
      siteInclusion = param('multi_clustering__reportree__site_inclusion')   
      extra = optional('multi_clustering__reportree__extra')
      nomenclature = !nomenclature_path.empty() ? "--nomenclature-file prev_nomenclature.tsv" : ""
      """
      reportree.py \
        -m ${metadata} \
        -align ${alignment} \
        -out 'gt' \
        --analysis grapetree \
        --columns_summary_report ${SUMMARY_COLUMNS} \
        --matrix-4-grapetree \
        --mx-transpose \
        --n_proc ${task.cpus} \
        --thr ${thr}  \
        --loci-called ${lociCalled} \
        --unzip \
        --site-inclusion ${siteInclusion} ${soi_par} ${zoom_par} ${zoom_subtree_par} ${extra} ${nomenclature} 
      touch gt_zooms.txt
      """
}

process reportree_hc {
    container "ghcr.io/genpat-it/reportree:2.4.1--088b6651b8"
    cpus { [32, params.max_cpus as int].min() }   
    input:
      path(alignment)
      path(metadata)
    output:
      path '*'
      path 'hc_*.nwk', emit: nwk_hc
      path 'hc_metadata_w_partitions.tsv', emit: metadata_hc
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'copy', "${params.outdir}/result", pattern: '{hc_*.tsv,*.txt,*.nwk}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "reportree_hc.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "reportree_hc.cfg" }
    script:
      hcMethod = param('multi_clustering__reportree__HC_threshold')
      lociCalled = param('multi_clustering__reportree__loci_called')   
      siteInclusion = param('multi_clustering__reportree__site_inclusion')   
      """
        reportree.py \
        -m ${metadata} \
        -align ${alignment} \
        -out hc \
        --analysis HC \
        --columns_summary_report ${SUMMARY_COLUMNS} \
        --mx-transpose \
        --loci-called ${lociCalled} \
        --site-inclusion ${siteInclusion} \
        --HC-threshold ${hcMethod}
      """
}

process find_closest {
    container "ghcr.io/genpat-it/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/multi_clustering__reportree:/scripts:ro"
    memory { taskMemory( 1.GB, task.attempt ) }
    input:
      path(matrix)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'sample_of_interest_summary.txt'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "find_closest.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "find_closest.cfg" }
    script:
      samples =  param('multi_clustering__reportree__sample_of_interest').split(',').collect { it.trim() }.join(',')
      threshold = param('multi_clustering__reportree__report_threshold')
      """
        /scripts/filter-matrix-distance.py ${matrix} ${threshold} ${samples} sample_of_interest_summary.txt
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
        cat ${metadata} | sed 's/${SAMPLE_COLUMN}/name/i' > augur_metadata.tsv
        METADATA_LIST=\$(head -n 1 augur_metadata.tsv | tr \$'\t' ' ')
        augur refine --tree ${nwk} --output-tree tree_tt.nwk --output-node-data refine.node.json --metadata augur_metadata.tsv
        augur export v2 --tree tree_tt.nwk --node-data refine.node.json --output auspice.json \
          --color-by-metadata \${METADATA_LIST} \
          --geo-resolutions ${GEO_RESOLUTION_COLUMNS} \
          --metadata augur_metadata.tsv \
          --lat-longs ${geodata}
      """
}


workflow multi_clustering__reportree {
    take: 
        input
        raw_metadata
        geodata
        nomenclature
    main:
        metadata = prepare_metadata(raw_metadata).metadata
        fastas = input.flatMap { it[1] }.collect()  
        alignment = maaft(fastas).alignment
        matrix = reportree_gt(alignment, metadata, nomenclature).matrix

        if (optional('multi_clustering__reportree__sample_of_interest') && optional('multi_clustering__reportree__report_threshold')) {
            find_closest(matrix)
        }
        if (optional('multi_clustering__reportree__HC_threshold')) {
          reportree_hc(alignment, metadata)
          augur(reportree_hc.out.nwk_hc, reportree_hc.out.metadata_hc, geodata)
        }
}

workflow {
  multi_clustering__reportree(getInput(),  param('metadata'), param('geodata'), optionalOrDefault('multi_clustering__reportree__nomenclature', getEmpty()));
}