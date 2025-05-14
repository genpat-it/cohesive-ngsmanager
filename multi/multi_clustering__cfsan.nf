nextflow.enable.dsl=2

include { taskMemory;flattenPath } from '../functions/common.nf'
include { getReferenceUnkeyed;getResult;getInput;param;optional;optionalBoolean } from '../functions/parameters.nf'

def ID_COLUMN = param('metadata_id_column')
def SUMMARY_DATE_ALIASES = param('metadata_date_aliases').tokenize(',\s').join('|')  
def GEO_RESOLUTION_COLUMNS = param('metadata_geo_column')  
def SUMMARY_COLUMNS = param('reportree__summary_columns')  

def METADATA_ID_FULL = param('metadata_with_full_id') as boolean

def IMAGES = [
  '2.2.1': 'staphb/cfsan-snp-pipeline:2.2.1',
  '2.0.2': 'cfsanbiostatistics/snp-pipeline@sha256:448787923371ade95217982814db25efb1e01287a8180b523d76a9f093f97d01'
]

def IQTREE_IMAGES = [
  'iqtree_2.3.6': 'quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0',
  'iqtree_1.6.12': 'quay.io/biocontainers/iqtree:1.6.12--he513fc3_0'
]

def DOCKER_IMAGE = IMAGES[param('cfsan__version')] ?: (exit 2, "params (cfsan__version) not valid");
def IQTREE_DOCKER_IMAGE = IQTREE_IMAGES[optional('mltree_method')] ?: ''

process cfsan_snp_pipeline {
    container DOCKER_IMAGE
    containerOptions = "-v ${workflow.projectDir}/scripts/multi_clustering__cfsan:/scripts:ro"
    input:
      path(summary)
      path(samples)
      tuple val(_), val(reference), path(refPath)
    output:
      path '**'
      path 'results/snpma.fasta', emit: snpma
      path 'results/snpma_preserved.fasta', emit: snpma_preserved
      path '{*.sh,*.log}', hidden: true
    stageInMode 'symlink'  
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'results/*.tsv', saveAs: { filename -> flattenPath(filename) }      
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'results/*.vcf', saveAs: { filename -> flattenPath(filename) }      
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'results/*.fasta', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "cfsan_snp_pipeline.cfg" }
    script:
      version = param('cfsan__version')
      extra_files_cleanup = optionalBoolean('cfsan__keep_intermediate_files') ? "" : "-o -name '*.bam' -o -name '*.pileup'"
      """
        trap "find  \\( -name '*.sam' ${extra_files_cleanup} \\) -delete" EXIT
        while IFS=\$'\t' read -r key file1 file2; do
            if (${METADATA_ID_FULL}); then
                KEY=\$key
            else
                KEY=\$(sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_.+/\\1/' <<< "\$file1")
            fi
            mkdir -p "samples/\${KEY}"
            [ -f "\$file1" ] && mv "\$file1" "samples/\${KEY}"/
            [ -f "\$file2" ] && mv "\$file2" "samples/\${KEY}"/
        done < ${summary}
        cfsan_snp_pipeline run -c /scripts/snppipeline_${version}.conf -m soft -o results --samples_dir samples ${refPath} >> cfsan_snp_pipeline.log
      """
}

process iqtree {
    container IQTREE_DOCKER_IMAGE
    stageInMode 'copy'
    input:
      path(snpma)
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
      if (optional('mltree_method') == 'iqtree_1.6.12') {
        """
          iqtree -s ${snpma} -nt AUTO
        """
      } else {    
        """
          iqtree \
            -s ${snpma} \
            -nt AUTO \
            -st DNA \
            -m MFP \
            --ufboot 1000 \
            --wbtl \
            --bnni
        """        
      }
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

workflow multi_clustering__cfsan {
    take: 
        input
        reference
        metadata
        geodata
    main:
        summary = input.collectFile { [ "summary.tsv", it[0] + '\t' + (it[1][0].getName() ?: '') + '\t' + (it[1][1]?.getName() ?: '') + '\n'] }
        files= input.collect { it[1] }

        cfsan_snp_pipeline(summary, files, reference)
        if (optionalBoolean('cfsan__filter_regions')) {
          snpma = cfsan_snp_pipeline.out.snpma_preserved
        } else {
          snpma = cfsan_snp_pipeline.out.snpma
        } 
        method = optional('mltree_method')
        if (method == 'fasttree') {
          nwk = fasttree(snpma).nwk
          augur(nwk, metadata, geodata, 'fasttree')
        } else if (IQTREE_DOCKER_IMAGE) {
          nwk = iqtree(snpma).nwk
          augur(nwk, metadata, geodata, 'iqtree')
        }
        if (optional('reportree__thr')) {
          reportree_gt(snpma, metadata)
        }
        if (optional('reportree__HC_method')) {
          nwk_hc = reportree_hc(snpma, metadata).nwk_hc
          augur_hc(nwk_hc, metadata, geodata, 'reportree')
        }
}

workflow {
    multi_clustering__cfsan(getInput(), getReferenceUnkeyed('fa'), param('metadata'), param('geodata'));
}