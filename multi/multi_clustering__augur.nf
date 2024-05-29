nextflow.enable.dsl=2

include { taskMemory } from '../functions/common.nf'
include { getInput;getReferenceUnkeyed;param;optional } from '../functions/parameters.nf'

def SUMMARY_DATE_ALIASES = param('multi_clustering__reportree__summary_date_aliases')  
def SAMPLE_COLUMN = param('multi_clustering__reportree__summary_sample_column')  
def GEO_RESOLUTION_COLUMNS = param('multi_clustering__reportree__summary_geo_column')  
def TRAITS_COLUMNS = param('multi_clustering__augur__traits_columns') 

process align {
    container "quay.io/biocontainers/augur:23.1.1--pyhdfd78af_1"
    memory { taskMemory( 5.GB, task.attempt ) }
    input:
      path(sequences)
      tuple val(_), val(ref_code), path(ref_path)
    output:
      path '*'
      path 'alignment.fasta', emit: alignment
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'alignment.fasta' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "align.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "align.cfg" }
    script:
      extra = optional('multi_clustering__augur__align_extra')  
      """
        for f in ${sequences}; do 
          cmp=`echo -n \$f | sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_[[:graph:]]+/\\1/'`
          awk "NR==1 {\\\$0=\\">\$cmp\\"}1" \$f > seq_\$f
        done     
        augur align \
        --nthreads auto \
        --sequences seq_* \
        --reference-sequence ${ref_path} \
        --remove-reference \
        --output alignment.fasta \
        --fill-gaps ${extra}  
      """
}

process tree {
    container "quay.io/biocontainers/augur:23.1.1--pyhdfd78af_1"
    input:
      path(alignment)
    output:
      path '*'
      path 'tree_raw.nwk', emit: tree_raw
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'tree_raw.nwk' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "tree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "tree.cfg" }
    script:
      method = param('multi_clustering__augur__tree_method')  
      extra = optional('multi_clustering__augur__tree_extra')  
      """
        augur tree \
          --alignment ${alignment} \
          --output tree_raw.nwk \
          --method ${method} \
          --nthreads auto ${extra}
      """
}

process refine {
    container "quay.io/biocontainers/augur:23.1.1--pyhdfd78af_1"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      path(tree)
      path(alignment)
      path(metadata)
    output:
      path '*'
      path 'tree.nwk', emit: tree
      path 'branch_lengths.json', emit: branch_lengths
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'tree.nwk' 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'branch_lengths.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "refine.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "refine.cfg" }
    script:
      coalescent = param('multi_clustering__augur__refine_coalescent')  
      extra = optional('multi_clustering__augur__refine_extra')  
      """
         augur refine \
            --tree ${tree} \
            --alignment ${alignment} \
            --metadata ${metadata} \
            --output-tree tree.nwk \
            --output-node-data branch_lengths.json \
            --timetree \
            --coalescent ${coalescent} \
            --date-confidence ${extra}
      """
}

process ancestral {
    container "quay.io/biocontainers/augur:23.1.1--pyhdfd78af_1"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      path(tree)
      path(alignment)
      tuple val(_), val(ref_code), path(ref_path)
    output:
      path '*'
      path 'nt_muts.json', emit: nt_muts
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'nt_muts.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "ancestral.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "ancestral.cfg" }
    script:
      inference = param('multi_clustering__augur__ancestral_inference')  
      extra = optional('multi_clustering__augur__ancestral_extra')      
      """
        augur ancestral \
            --tree ${tree} \
            --alignment ${alignment} \
            --root-sequence ${ref_path} \
            --output-node-data nt_muts.json \
            --inference ${inference} ${extra}
      """
}

process translate {
    container "quay.io/biocontainers/augur:23.1.1--pyhdfd78af_1"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      path(tree)
      path(nt_muts)
      tuple val(_), val(ref_code), path(ref_path)
    output:
      path '*'
      path 'aa_muts.json', emit: aa_muts
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'aa_muts.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "translate.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "translate.cfg" }
    script:
      extra = optional('multi_clustering__augur__translate_extra')      
      """
        augur translate \
            --tree ${tree} \
            --ancestral-sequences ${nt_muts} \
            --reference-sequence ${ref_path} \
            --output aa_muts.json ${extra}
      """
}

process traits {
    container "quay.io/biocontainers/augur:23.1.1--pyhdfd78af_1"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      path(tree)
      path(metadata)
    output:
      path '*'
      path 'traits.json', emit: traits
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'traits.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "traits.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "traits.cfg" }
    script:
      columns = TRAITS_COLUMNS   
      extra = optional('multi_clustering__augur__traits_extra')      
      """
        augur traits \
            --tree ${tree} \
            --metadata ${metadata} \
            --output-node-data traits.json \
            --columns ${columns} \
            --confidence ${extra}
      """
}

process prepare_metadata {
    container "quay.io/biocontainers/augur:23.1.1--pyhdfd78af_1"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      path(metadata)
    output:
      path '*'
      path 'augur_metadata.tsv', emit: metadata
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'augur_metadata.tsv' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "prepare_metadata.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "prepare_metadata.cfg" }
    script:
      """
         cat ${metadata} | sed 's/${SAMPLE_COLUMN}/name/i' | sed -E 's/${SUMMARY_DATE_ALIASES}/date/i' > augur_metadata.tsv 
      """
}

process export {
    container "quay.io/biocontainers/augur:23.1.1--pyhdfd78af_1"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      path(tree)
      path(metadata)
      path(branch_lengths)
      path(traits)
      path(nt_muts)
      path(aa_muts)
      path(lat_longs)
    output:
      path '*'
      path 'auspice.json', emit: auspice
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'auspice.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "export.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "export.cfg" }
    script:
      extra = optional('multi_clustering__augur__export_extra')      
      """
        METADATA_LIST=\$(head -n 1 ${metadata} | tr \$'\t' ' ')
        augur export v2 \
            --tree ${tree} \
            --metadata ${metadata} \
            --node-data ${branch_lengths} ${traits} ${nt_muts} ${aa_muts} \
            --lat-longs ${lat_longs} \
            --color-by-metadata \${METADATA_LIST} \
            --geo-resolutions ${GEO_RESOLUTION_COLUMNS} \
            --output auspice.json ${extra}
      """
}

workflow multi_clustering__augur {
    take: 
      reference
      raw_metadata
      geodata
      ref2
    main:  
      metadata = prepare_metadata(raw_metadata).metadata
      fastas = getInput().flatMap { it[1] }.collect()
      alignment = align(fastas, reference).alignment
      tree_raw = tree(alignment).tree_raw
      tree = refine(tree_raw, alignment, metadata).tree
      nt_muts = ancestral(tree, alignment, reference).nt_muts
      aa_muts = translate(tree, nt_muts, reference).aa_muts
      traits = traits(tree, metadata).traits
      export(tree, metadata, refine.out.branch_lengths, traits, nt_muts, aa_muts, geodata)
  }

workflow {
    reference = getReferenceUnkeyed('gb')      
    multi_clustering__augur(reference, param('metadata'), param('geodata'), reference)
}
