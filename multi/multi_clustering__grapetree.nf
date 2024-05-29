nextflow.enable.dsl=2

include { parseRISCD;taskMemory } from '../functions/common.nf'
include { _getAlleles;param } from '../functions/parameters.nf'

def GEO_RESOLUTION_COLUMNS = param('multi_clustering__reportree__summary_geo_column')
def SUMMARY_DATE_ALIASES = param('multi_clustering__reportree__summary_date_aliases')
def SAMPLE_COLUMN = param('multi_clustering__reportree__summary_sample_column')

process extract_cgMLST {
    container "${LOCAL_REGISTRY}/bioinfo/chewbbaca-w-schemas:2.8.5--7742d1fae0"
    memory { taskMemory( 5.GB, task.attempt ) }
    input:
      path(alleles)
    output:
      path '**'
      path 'cgMLST.tsv', emit: cgMLST
      path '*.sh', hidden: true
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'mdata_stats.tsv', saveAs: { "missing_loci.tsv" }
    publishDir mode: 'copy', "${params.outdir}/result", pattern: 'cgMLST.tsv', saveAs: { "cgMLST.tsv" }
    publishDir mode: 'copy', "${params.outdir}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "extract_cgMLST.cfg" }
    script:
      """
        for file in ${alleles} ; do awk 'FNR==1{print ""}1' \${file} | sed 's/,/\t/g' | sed -E "s/^[^SF][^ai]\\S+/\${file}/" | sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_[[:graph:]]+/\\1/'; done | sort -ru  > results_alleles_all.tsv
        chewie ExtractCgMLST -i results_alleles_all.tsv -o . > extract_cgMLST.log
      """
}

process dists {
    container "quay.io/biocontainers/cgmlst-dists:0.4.0--hec16e2b_2"
    memory { taskMemory( 5.GB, task.attempt ) }
    input:
      path(cgMLST)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{*.csv}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "cgmlst-dists.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "cgmlst-dists.cfg" }
    script:
      """
        cgmlst-dists -c ${cgMLST} > cgMLST_dists_matrix.csv
      """
}

process grapetree {
    container "quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0"
    memory { taskMemory( 5.GB, task.attempt ) }
    input:
      path(cgMLST)
    output:
      path '*'
      path '*.sh', hidden: true
      path("cgMLST_NJ.nwk"), emit: nwk_nj
    publishDir mode: 'copy', "${params.outdir}/result", pattern: '*.nwk'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "grapetree.cfg" }
    script:
      """
        grapetree -p ${cgMLST} > cgMLST.nwk
        grapetree --method RapidNJ -p ${cgMLST} > cgMLST_NJ.nwk
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
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log', saveAs: { "augur.log" }
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

workflow multi_clustering__grapetree {
    take: 
        input
        metadata
        geodata
    main:
        cgMLST = extract_cgMLST(input.collect()).cgMLST
        dists(cgMLST)
        nwk_nj = grapetree(cgMLST).nwk_nj
        augur(nwk_nj, metadata, geodata)
}

workflow {
    multi_clustering__grapetree(getInput(), param('metadata'), param('geodata'));
}

def getInput() {
    if (!params.containsKey('input')) {
      exit 2, "missing required param: input";
    }
    def schema = params.containsKey('schema') ? params.schema : null
    assert params.input instanceof ArrayList
    params.input.inject(Channel.empty()) {
        res, val -> res.mix(_getAlleles(val.cmp, val.riscd, schema))
    }        
}

