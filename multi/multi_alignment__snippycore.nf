nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { param;optionalBoolean;optional;getInputOf;getReferenceUnkeyed;getInputFolders;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent  } from '../functions/parameters.nf'

def ex = executionMetadata()

def ENTRYPOINT = "multi_alignment__snippycore"

def ID_COLUMN = param('metadata_id_column')
def SUMMARY_DATE_ALIASES = param('multi_clustering__reportree__summary_date_aliases')  
def GEO_RESOLUTION_COLUMNS = param('multi_clustering__reportree__summary_geo_column')  
def SUMMARY_COLUMNS = param('multi_clustering__reportree__summary_columns')  


process snippy {
    container "ghcr.io/genpat-it/snippy:4.5.1--7be4a1c45a"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 8.GB, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(reference_path)
    output:
      path 'snippy', emit: results
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base_ref}.log" }    
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_snippy_${reference}"
      is_fasta = r1.getName() ==~ /.+\.fa(sta)?$/
      if (is_fasta) {
        """
        trap "find -name "*.?am" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --ctgs ${r1} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        """     
      } else if (isIlluminaPaired(reads)) {
        """
        trap "find -name "*.?am" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --R1 ${r1} --R2 ${r2} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        """
      } else if (isIonTorrent(reads)) {
        """
        trap "find -name "*.?am" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --se ${r1} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        """      
      }      
}

process snippy_core {
    container "ghcr.io/genpat-it/snippy:4.5.1--7be4a1c45a"
    memory { taskMemory( 4.GB, task.attempt ) }
    input:
      path vcf_folders, stageAs: 'data?'
      tuple val(ref_riscd), val(ref_code), path(ref_file)
    output:
      path '*'
      path 'core.aln', emit: snpma
      path 'clean.full.aln', emit: full_alignment
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{core*}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "snippy_core.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log', saveAs: { "snippy_core.log" }
    script:
      // XXX renaming folders getting sample name from the first vcf file inside
      """          
        #!/bin/bash -euo pipefail
        mkdir inputs && cd inputs && for dir in ${vcf_folders}; do ln -s ../\${dir} `ls ../\${dir}/*.vcf | head -n 1 | sed -E 's/.+DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_[[:graph:]]+/\\1/'`; done && cd ..
        snippy-core --ref ${ref_file} --prefix core --inprefix snps inputs/*
        snippy-clean_full_aln core.full.aln > clean.full.aln
      """
}

process iqtree {
    container 'quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0'
    stageInMode 'copy'
    input:
      path(snpma)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path '*.treefile', emit: nwk
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.ufboot', saveAs: { "iqtree.ufboot" } 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.mldist', saveAs: { "iqtree.mldist" } 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.treefile', saveAs: { "iqtree.nwk" } 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.iqtree', saveAs: { "iqtree_report.txt" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "iqtree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "iqtree.cfg" }
    script:
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

process fasttree {
    container 'quay.io/biocontainers/fasttree:2.1.11--h031d066_4'
    input:
      path(snpma)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path '*.nwk', emit: nwk
    publishDir mode: 'copy', "${params.outdir}/result", pattern: '*.nwk'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "fasttree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "fasttree.cfg" }
    script:
      """
        FastTree -nt -gtr -gamma  ${snpma} > snpma.nwk
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

process gubbins {  
    container "quay.io/biocontainers/gubbins:2.2.1--py35_0"
    input:
      path(alignment)
    output:
      path '*'
      path 'gubbins.filtered_polymorphic_sites.fasta', emit: clean_alignment
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "gubbins.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "gubbins.cfg" }
    script:
      """
      run_gubbins.py -p gubbins ${alignment}
      """
}

process reportree_gt {
    container "ghcr.io/genpat-it/reportree:2.4.1--088b6651b8"
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
      thr = param('multi_clustering__reportree__thr')      
      siteInclusion = param('multi_clustering__reportree__site_inclusion')   
      lociCalled = param('multi_clustering__reportree__loci_called')   
      extra = optional('multi_clustering__reportree__extra')
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
        --site-inclusion ${siteInclusion} ${extra}
      touch gt_zooms.txt
      """
}

workflow multi_alignment__snippycore {
    take: 
        reads
        pre_calculated_result_folders      
        reference
        metadata
        geodata        
    main:
        reads.combine(reference)
                .multiMap { 
                    reads: it[0..1] // riscd, R[]
                    reference:  it[2..4] // riscd, code, path
                }.set { input }
        folders = snippy(input.reads,input.reference).results

        vcf_folders = folders.mix(pre_calculated_result_folders.map { it[1] }).collect()
        
        if (param('multi_alignment__snippycore__alignment_type') == 'core') {
          ma_raw = snippy_core(vcf_folders, reference).snpma      
        } else {
          ma_raw = snippy_core(vcf_folders, reference).full_alignment
        }

        if (optionalBoolean('multi_alignment__snippycore__remove_recombinations')) {
          ma = gubbins(ma_raw).clean_alignment
        } else {
          ma = ma_raw
        } 

        method = param('multi_alignment__snippycore__tree_method')
        if (method == 'fasttree') {
          nwk = fasttree(ma).nwk
          augur(nwk, metadata, geodata)
        } else if (method == 'iqtree') {
          nwk = iqtree(ma).nwk
          augur(nwk, metadata, geodata)
        }
        reportree_gt(ma, metadata)
}

workflow {
    // use getInputOf instead of getInput to avoid warnings
    reads = getInputOf('{*.vcf,*.fastq*}').filter( ~/^.*\/1PP_.*/ )
    pre_calculated_result_folders = getInputFolders().filter( ~/^.*\/2AS_.*/ )

    // pre_calculated_result_folders have to be calculated using the same reference and format (gb or fasta)
    multi_alignment__snippycore(reads, pre_calculated_result_folders, getReferenceUnkeyed('any'), param('metadata'), param('geodata'))
}