nextflow.enable.dsl=2

include { param;optional;optWrap;optionalBoolean } from '../functions/parameters.nf'
include { taskMemory;taskCpus } from '../functions/common.nf'

process iqtree {
    container 'quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0'
    memory { taskMemory( 50.GB, task.attempt ) }
    stageInMode 'copy'
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.ufboot', saveAs: { "iqtree.ufboot" } 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.mldist', saveAs: { "iqtree.mldist" } 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.treefile', saveAs: { "iqtree.nwk" } 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.iqtree', saveAs: { "iqtree_report.txt" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "iqtree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "iqtree.cfg" }
    input:
      path(snpma)
      val(polymorphic_sites_only)
      path(constant_sites)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path '*.treefile', emit: nwk
    script:
      if (polymorphic_sites_only && constant_sites.empty()) {       
        model = 'MFP+ASC'         
      } else {
        model = 'MFP'
      }
      // fconst = ''
      fconst = !(constant_sites.empty()) ? "-fconst `cat ${constant_sites}`" : ''
      pathogen = optionalBoolean('iqtree__pathogen') ?  '--pathogen' : ''
      extra = optional('mltree_extra')
      """      
        iqtree ${fconst} \
          -s ${snpma} \
          -nt AUTO \
          -st DNA \
          -m ${model} \
          --ufboot 1000 \
          --wbtl \
          --bnni ${pathogen} ${extra}
      """
}

process augur {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 4.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'auspice.json', saveAs: { "auspice_${source}.json" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "augur_${source}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "augur_${source}.cfg" }
    input:
      path(nwk)
      path(metadata)
      path(geodata)
      val(source)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    script:
        ID_COLUMN = param('metadata_id_column')
        SUMMARY_DATE_ALIASES = param('metadata_date_aliases').tokenize(',\s').join('|')
        GEO_RESOLUTION_COLUMNS = param('metadata_geo_column')          
      """
        sed 's/${ID_COLUMN}/name/i' ${metadata} | sed -E 's/${SUMMARY_DATE_ALIASES}/date/i' > augur_metadata.tsv
        METADATA_LIST=\$(head -n 1 augur_metadata.tsv | tr \$'\t' ' ')
        augur refine --tree ${nwk} --output-tree tree_tt.nwk --output-node-data refine.node.json --metadata augur_metadata.tsv
        augur export v2 --tree tree_tt.nwk --node-data refine.node.json --output auspice.json \
          --color-by-metadata \${METADATA_LIST} \
          --geo-resolutions ${GEO_RESOLUTION_COLUMNS} \
          --metadata augur_metadata.tsv \
          --lat-longs ${geodata}
      """
}

process fasttree {
    container 'quay.io/biocontainers/fasttree:2.1.11--h031d066_4'
    memory { taskMemory( 10.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.nwk', saveAs: { "fasttree.nwk" }      
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "fasttree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "fasttree.cfg" }
    input:
      path(snpma)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path '*.nwk', emit: nwk
    script:
      extra = optional('mltree_extra')
      """
        FastTree -nt -gtr -gamma ${extra} ${snpma} > snpma.nwk
      """
}

process reportree_alignment {
    container "${params.local_registry}/bioinfo/reportree:2.5.4--0f5f86c689"
    memory { taskMemory( 150.GB, task.attempt ) }
    cpus { taskCpus(96, task.attempt) } 
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'gt_zooms.txt', saveAs: { "zooms.txt" }
    publishDir mode: 'copy', "${params.outdir}/result", pattern: '{*.tsv,*.txt,*.nwk}'
    publishDir mode: 'copy', "${params.outdir}/result", pattern: '*_*[0-9]'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "reportree_${prefix}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "reportree_${prefix}.cfg" }
    input:
      path(alignment)
      path(metadata)
      val(analysis)
      path(nomenclature_path), stageAs: 'prev_nomenclature.tsv'
    output:
      path '*'
      path "${prefix}_dist_hamming.tsv", emit: matrix
      path "${prefix}*.nwk", emit: nwk
      path "${prefix}_zooms.txt", emit: zooms, optional: true
      path "${prefix}_metadata_w_partitions.tsv", emit: metadata_w_partitions
      path "${prefix}_SAMPLES_OF_INTEREST_partitions_summary.tsv", emit: soi_partition_summary, optional: true
      path '{*.sh,*.log}', hidden: true
    script:
      soi = optional('reportree__sample_of_interest').replaceAll(/[\s\n\t\r"'$\{]/, "")	
      soi_par = optWrap('reportree__sample_of_interest', "--sample_of_interest ${soi}")
      zoom_par = optWrap('reportree__zoom_cluster_of_interest', '--zoom-cluster-of-interest {}')
      zoom_subtree_par = optWrap('reportree__subtree_of_interest', '--subtree-of-interest {}')
      site_inclusion = optWrap('reportree__site_inclusion', '--site-inclusion {}') 
      sample_ATCG_content = optWrap('reportree__sample_atcg_content', '--sample-ATCG-content {}') 
      extra = optional("reportree__${analysis}_extra")
      nomenclature = nomenclature_path.exists() && !nomenclature_path.empty() ? "--nomenclature-file prev_nomenclature.tsv" : ""
      if (analysis == 'grapetree') {
        prefix = 'gt'
        thr = optWrap('reportree__thr', '--thr {}')      
        nproc = "--n_proc ${task.cpus}"
        matrix_grapetre = '--matrix-4-grapetree'
        hc_method = ''
        hc_threshold = ''
        pct_hc_threshold = ''
      } else {
        prefix = 'hc'
        thr = ''
        nproc = ''
        matrix_grapetre = ''
        // hc_method,hc_threshold,pct_hc_threshold mutually exclusive
        // just one method at once
        hc_method = optWrap('reportree__HC_method', '--HC-threshold {}')
        hc_threshold = optWrap('reportree__HC_threshold', '--HC-threshold {}')
        pct_hc_threshold = optWrap('reportree__pct_HC_threshold', '--pct-HC-threshold {}')
      }
      SUMMARY_COLUMNS = param('reportree__summary_columns')
      """
      sed 's/-/N/g' ${alignment} > preprocessed_alignment.fasta
      reportree.py \
        -m ${metadata} \
        -align preprocessed_alignment.fasta \
        -out '${prefix}' \
        --analysis ${analysis} \
        --columns_summary_report ${SUMMARY_COLUMNS} \
        ${nproc} \
        ${thr}  \
        ${sample_ATCG_content} \
        ${site_inclusion} \
        --unzip \
        ${hc_method} ${hc_threshold} ${pct_hc_threshold} \
        ${matrix_grapetre} ${soi_par} ${zoom_par} ${zoom_subtree_par} ${nomenclature} ${extra} 
      """
}

process prepare_reportree_metadata {
    container "ubuntu:20.04"
    memory { taskMemory( 2.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'reportree_metadata.tsv' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "prepare_metadata.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "prepare_metadata.cfg" }
    input:
      path(metadata)
    output:
      path '*'
      path 'reportree_metadata.tsv', emit: metadata
    script:
      soi = optional('reportree__sample_of_interest').replaceAll(/[\s\n\t\r"'$\{]/, "")	
      SUMMARY_DATE_ALIASES = param('metadata_date_aliases').tokenize(',\s').join('|') 
      """
        awk -v SOI="${soi}," '
          BEGIN { FS = OFS = "\\t" }
          {
              if (NR == 1) {
                  # Add a new "group" column to the header
                  \$1 = \$1 FS "group"
              } else {
                  # Determine the group based on conditions
                  if (index(SOI, \$1 ","))
                      group = "sample of interest"
                  else if (\$1 ~ /20[0-9][0-9]\\.EXT\\.[0-9].*/)
                      group = "imported"
                  else if (\$1 ~ /20[0-9][0-9]\\.[A-Z]+\\.[0-9].*/)
                      group = "internal"
                  else
                      group = "external"
                  # Append the group to the current line
                  \$1 = \$1 FS group
              }
          } 1
        ' ${metadata} | sed -E 's/${SUMMARY_DATE_ALIASES}/date/i' > reportree_metadata.tsv	
      """     
}

process find_closest {
    container "${params.local_registry}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions "-v ${workflow.projectDir}/scripts/multi_clustering__reportree:/scripts:ro"
    memory { taskMemory( 1.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'sample_of_interest_summary.txt'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "find_closest.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "find_closest.cfg" }
    input:
      path(matrix)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    script:
      samples =  param('reportree__sample_of_interest').split(',').collect { it.trim() }.join(',')
      threshold = param('closest_samples_report_threshold')
      """
        /scripts/filter-matrix-distance.py ${matrix} ${threshold} ${samples} sample_of_interest_summary.txt
      """
}