nextflow.enable.dsl=2

include { getInput;param;optional } from '../functions/parameters.nf'
include { taskMemory;flattenPath } from '../functions/common.nf'

def ID_COLUMN = param('metadata_id_column')
def SUMMARY_DATE_ALIASES = param('metadata_date_aliases').tokenize(',\s').join('|')  
def GEO_RESOLUTION_COLUMNS = param('metadata_geo_column')  
def SUMMARY_COLUMNS = param('reportree__summary_columns')  

def METADATA_ID_FULL = param('metadata_with_full_id') as boolean

def IMAGES = [
  '4.1': 'quay.io/staphb/ksnp4:4.1',
  '3': "${LOCAL_REGISTRY}/bioinfo/ksnp3:3.0--addd2c2d0e"
]
def KSNP_EXECUTABLES = [
  '4.1': 'kSNP4',
  '3': 'kSNP3'
]

def KSNP_DOCKER_IMAGE = IMAGES[param('ksnp__version')] ?: (exit 2, "params (ksnp__version) not valid");
def KSNP_EXECUTABLE = KSNP_EXECUTABLES[param('ksnp__version')] ?: (exit 2, "params (ksnp__version) not valid");

process kchooser {
    container 'quay.io/staphb/ksnp4:4.1'
    containerOptions = "-u 0:0"
    input:
      path(summary)
      path(assembly)
    output:
      path('kmers_size.txt'), emit: ksize
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'copy', "${params.outdir}/result", pattern: "Kchooser*_input.report"  
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "kchooser.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "kchooser.cfg" }
    script:
      """
        while IFS=\$'\t' read -r key file; do
            if (${METADATA_ID_FULL}); then
                KEY=\$key
            else
                KEY=\$(sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([[:digit:]]+\\.[[:alpha:]]+\\.[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+).+/\\1/' <<< "\$file" | tr '.' '-')
            fi
            mv \$file \${KEY}.fasta
            echo -e "\$(pwd)/\${KEY}.fasta\t\$KEY" 
        done < ${summary} > input.tsv
        Kchooser4 -in input.tsv
        grep -oP 'optimum value of k is \\K\\d+' Kchooser*_input.report > kmers_size.txt
      """
}

process ksnp {
    container KSNP_DOCKER_IMAGE
    containerOptions = "-u 0:0"
    input:
      path(summary)
      path(assembly)
      val(kmers_size)
      val(analysis_type)
    output:
      path("${outfile}"), emit: fasta
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'copy', "${params.outdir}/result", pattern: "${outfile}", saveAs: { filename -> flattenPath(filename) }   
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "ksnp.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "ksnp.cfg" }
    script:
      extra_params = (analysis_type == 'core' ? '-core' : '')
      outfile = (analysis_type == 'core' ? 'core_SNPs_matrix.fasta' : 'SNPs_all_matrix.fasta')
      """
        trap "rm -Rf results/TemporaryFilesToDelete" EXIT
        while IFS=\$'\t' read -r key file; do
            if (${METADATA_ID_FULL}); then
                KEY=\$key
            else
                KEY=\$(sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([[:digit:]]+\\.[[:alpha:]]+\\.[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+).+/\\1/' <<< "\$file" | tr '.' '-')
            fi
            mv \$file \${KEY}.fasta
            echo -e "\$(pwd)/\${KEY}.fasta\t\$KEY" 
        done < ${summary} > input.tsv       
        ${KSNP_EXECUTABLE} -in input.tsv -k ${kmers_size} ${extra_params} -outdir results
        mv results/${outfile} results/${outfile}.tmp && awk '/^>/ {gsub("-", ".")} {print}' results/${outfile}.tmp > ${outfile}
      """
}

process iqtree {
    container 'quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0'
    stageInMode 'copy'
    input:
      path(fasta)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path "*.treefile", emit: nwk
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.ufboot', saveAs: { "iqtree.ufboot" } 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.mldist', saveAs: { "iqtree.mldist" } 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.treefile', saveAs: { "iqtree.nwk" } 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.iqtree', saveAs: { "iqtree_report.txt" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "iqtree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "iqtree.cfg" }
    script:
      """
        iqtree \
          -s ${fasta} \
          -nt AUTO \
          -st DNA \
          -m MFP \
          --ufboot 1000 \
          --wbtl \
          --bnni
        sed -i 's/-/./g' *.treefile
      """
}

process augur {
    container "quay.io/biocontainers/augur:22.0.0--pyhdfd78af_0"
    memory { taskMemory( 4.GB, task.attempt ) }
    input:
      path(nwk)
      path(metadata)
      path(geodata)
      val(source)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'auspice.json', saveAs: { "auspice_${source}.json" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "augur_${source}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "augur_${source}.cfg" }
    script:
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
    input:
      path(snpma)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path '*.nwk', emit: nwk
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.nwk', saveAs: { "fasttree.nwk" }      
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "fasttree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "fasttree.cfg" }
    script:
      """
        FastTree -nt -gtr -gamma  ${snpma} > snpma.nwk
      """
}

process reportree_gt {
    container "${LOCAL_REGISTRY}/bioinfo/reportree:2.5.4--0f5f86c689"
    cpus 64     
    input:
      path(alignment)
      path(metadata)
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
      thr = param('reportree__thr')      
      soi = optional('reportree__sample_of_interest')
      soi_par = soi ? " --sample_of_interest ${soi}" : ''
      zoom = optional('reportree__zoom_cluster_of_interest')
      zoom_par = zoom ? " --zoom-cluster-of-interest  ${zoom}" : ''
      zoom_subtree = optional('reportree__subtree_of_interest')
      zoom_subtree_par = zoom_subtree ? " --subtree-of-interest  ${zoom_subtree}" : ''
      lociCalled = param('reportree__loci_called')   
      siteInclusion = param('reportree__site_inclusion')   
      extra = optional('reportree__extra')
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
        --site-inclusion ${siteInclusion} ${soi_par} ${zoom_par} ${zoom_subtree_par} ${extra}
      touch gt_zooms.txt
      """
}

process reportree_hc {
    container "${LOCAL_REGISTRY}/bioinfo/reportree:2.5.4--0f5f86c689"
    cpus 64     
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
      hcMethod = param('reportree__HC_method')
      lociCalled = param('reportree__loci_called')   
      siteInclusion = param('reportree__site_inclusion')   
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

workflow augur_hc {
    take: 
        nwk
        metadata
        geodata
        source
    main:
      augur(nwk,metadata,geodata,source)
}

workflow multi_clustering__ksnp {
    take: 
        input
        kmers_size
        analysis_type
        metadata
        geodata
    main:
        summary = input.collectFile { [ "summary.tsv", it[0] + '\t' + it[1][0].getName() + '\n' ] }
        files= input.collect { it[1] }

        if (!kmers_size) {
          kmers_size = kchooser(summary, files).ksize.first().map { it.text.trim() }
        }

        ma = ksnp(summary, files, kmers_size, analysis_type).fasta        
        method = optional('mltree_method')
        if (method == 'fasttree') {
          nwk = fasttree(ma).nwk
          augur(nwk, metadata, geodata, 'fasttree')
        } else if (method == 'iqtree') {
          nwk = iqtree(ma).nwk
          augur(nwk, metadata, geodata, 'iqtree')
        }
        if (optional('reportree__thr')) {
          reportree_gt(ma, metadata)
        }
        if (optional('reportree__HC_method')) {
          nwk_hc = reportree_hc(ma, metadata).nwk_hc
          augur_hc(nwk_hc, metadata, geodata, 'reportree')
        }
}

workflow {
    multi_clustering__ksnp(getInput(), optional('kmers_size'), param('analysis_type'), param('metadata'), param('geodata'));
}
