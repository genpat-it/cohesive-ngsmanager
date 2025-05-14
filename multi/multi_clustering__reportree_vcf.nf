nextflow.enable.dsl=2

include { taskMemory;getEmpty;filename } from '../functions/common.nf'
include { _getAlleles;param;optional;optionalBoolean;optionalOrDefault;getVCFs } from '../functions/parameters.nf'

def SUMMARY_DATE_ALIASES = param('multi_clustering__reportree__summary_date_aliases').tokenize(',\s').join('|')  
def SUMMARY_COLUMNS = param('multi_clustering__reportree__summary_columns')  
def ID_COLUMN = param('metadata_id_column')  
def GEO_RESOLUTION_COLUMNS = param('multi_clustering__reportree__summary_geo_column')  

def METADATA_ID_FULL = param('metadata_with_full_id') as boolean

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
      soi = optional('multi_clustering__reportree__sample_of_interest').replaceAll(/[\s\n\t\r"'$\{]/, "")	
      if (soi) {
        """
          awk -v SOI="${soi}," 'BEGIN{ FS=OFS="\\t" } {\$1 = \$1 FS (NR==1? "group" : (index(SOI, \$1",")? "sample of interest" : "other")) }1' ${metadata} | sed -E 's/${SUMMARY_DATE_ALIASES}/date/i' > reportree_metadata.tsv	
        """
      } else {
        """
          sed -E 's/${SUMMARY_DATE_ALIASES}/date/i' ${metadata} > reportree_metadata.tsv
        """
      }
}

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

process reportree_gt {
    container "${LOCAL_REGISTRY}/bioinfo/reportree:2.5.4--0f5f86c689"
    cpus 64     
    input:
      path(summary)
      path(vcfs)
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
      extra = optional('multi_clustering__reportree__extra')
      nomenclature = !nomenclature_path.empty() ? "--nomenclature-file prev_nomenclature.tsv" : ""
      """
      #!/bin/bash -euo pipefail

      while IFS=\$'\t' read -r key file; do
          if (${METADATA_ID_FULL}); then
              KEY=\$key
          else
              KEY=\$(sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_.+/\\1/' <<< "\$file")
          fi
          mv \$file \${KEY}
          echo "\${KEY}" 
      done < ${summary} > vcf_files.txt
     
      reportree.py \
        -m ${metadata} \
        -vcf vcf_files.txt \
        -out 'gt' \
        --analysis grapetree \
        --columns_summary_report ${SUMMARY_COLUMNS} \
        --matrix-4-grapetree \
        --mx-transpose \
        --n_proc ${task.cpus} \
        --thr ${thr}  \
        --unzip ${soi_par} ${zoom_par} ${zoom_subtree_par} ${extra} ${nomenclature} 
      touch gt_zooms.txt        
      """      
}

process reportree_hc {
    container "${LOCAL_REGISTRY}/bioinfo/reportree:2.5.4--0f5f86c689"
    cpus 64     
    input:
      path(summary)
      path(vcfs)
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

      reportree.py \
        -m ${metadata} \
        -vcf vcf_files.txt \
        -out hc \
        --analysis HC \
        --columns_summary_report ${SUMMARY_COLUMNS} \
        --mx-transpose \
        --HC-threshold ${hcMethod}
      """
}

process find_closest {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
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
        cat ${metadata} | sed 's/${ID_COLUMN}/name/i' | sed -E 's/${SUMMARY_DATE_ALIASES}/date/i' > augur_metadata.tsv
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

        if (optionalBoolean('multi_clustering__reportree__vcf_normalization')) {
          vcf_input = bcftools(input).filtered_vcf
        } else {
          vcf_input = input
        }

        summary = vcf_input.collectFile { [ "summary.tsv", it[0] + '\t' + filename(it[1]) + '\n' ] }
        files= vcf_input.collect { it[1] }

        matrix = reportree_gt(summary, files, metadata, nomenclature).matrix

        if (optional('multi_clustering__reportree__sample_of_interest') && optional('multi_clustering__reportree__report_threshold')) {
            find_closest(matrix)
        }
        if (optional('multi_clustering__reportree__HC_threshold')) {
          reportree_hc(summary, files, metadata)
          augur(reportree_hc.out.nwk_hc, reportree_hc.out.metadata_hc, geodata)
        }
}

workflow {
  multi_clustering__reportree(getVCFs(),  param('metadata'), param('geodata'), optionalOrDefault('multi_clustering__reportree__nomenclature', getEmpty()));
}