nextflow.enable.dsl=2

include { taskMemory;flattenPath } from '../functions/common.nf'
include { getReferenceUnkeyed;getResult;getInput;param;optional;optionalBoolean } from '../functions/parameters.nf'

def ID_COLUMN = param('metadata_id_column')
def SUMMARY_DATE_ALIASES = param('metadata_date_aliases').tokenize(',\s').join('|')  
def GEO_RESOLUTION_COLUMNS = param('metadata_geo_column')  
def SUMMARY_COLUMNS = param('reportree__summary_columns')  

def METADATA_ID_FULL = param('metadata_with_full_id') as boolean

process mashtree {
    container 'quay.io/biocontainers/mashtree:1.4.6--pl5321h031d066_0'
    cpus params.max_cpus
    input:
      path(summary)
      path(samples)
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
      path 'mashtree.nwk', emit: nwk 
      path 'mashtree_distance.tsv', emit: matrix_distance 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'mashtree.nwk'
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'mashtree_distance.tsv'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.command.out', saveAs: { "mashtree.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "mashtree.cfg" }
    script:
      accurate_params = optionalBoolean('mashtree__accurate') ? '--mindepth 0' : ''
      extra = optional('mashtree__extra')
      """
        while IFS=\$'\t' read -r key file1 file2; do
            if (${METADATA_ID_FULL}); then
                KEY=\$key
            else
                KEY=\$(sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([[:digit:]]+\\.[[:alpha:]]+\\.[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+).+/\\1/' <<< "\$file1")
            fi
            if [ -f "\$file1" ] && [[ "\$file1" == *.f*a ]]; then
               mv "\$file1" "\${KEY}.fasta" && echo -e "\${KEY}.fasta"      
            fi       
        done < ${summary} > input_list.txt
        mashtree --numcpus ${task.cpus} --file-of-files input_list.txt --outmatrix mashtree_distance.tsv --outtree mashtree.nwk ${accurate_params} ${extra} 
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

workflow multi_clustering__mashtree {
    take: 
        input
        metadata
        geodata
    main:
        summary = input.collectFile { [ "summary.tsv", it[0] + '\t' + (it[1][0].getName() ?: '') + '\t' + (it[1][1]?.getName() ?: '') + '\n'] }
        files= input.collect { it[1] }

        matrix_distance = mashtree(summary, files).matrix_distance
        augur(mashtree.out.nwk, metadata, geodata, 'mashtree')
}

workflow {
    multi_clustering__mashtree(getInput(), param('metadata'), param('geodata'));
}