include { taskMemory } from '../functions/common.nf'
include { optWrap;optionalOrDefault;optional } from '../functions/parameters.nf'

process augur_filter {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 4.GB, task.attempt ) }
    cpus 2
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{*.fasta,*.tsv}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "augur_filter.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "augur_filter.cfg" }
    input:
      path(sequences)
      path(seq_index)
      path(metadata)
    output:
      path 'filtered_metadata.tsv', emit: filtered_metadata
      path 'filtered_sequences.fasta', emit: filtered_sequences
      path '{*.sh,*.log}', hidden: true
    script:
      query = optWrap('_augur__filter_query', '--query "{}"')    
      extra = optional('augur__filter_extra')
      min_length = optWrap('augur__filter_min_length', '--min-length "{}"')    
      max_length = optWrap('augur__filter_max_length', '--max-length "{}"')   
      metadata_id_column = optionalOrDefault('metadata-id-columns', "Id name strain") 
      """
        augur filter \
            --sequences ${sequences} \
            --metadata ${metadata} \
            --sequence-index ${seq_index} \
            --metadata-id-columns ${metadata_id_column} \
            ${min_length} ${max_length} ${query} ${extra} \
            --output filtered_sequences.fasta \
            --output-metadata filtered_metadata.tsv 
      """
}

process augur_filter_subsampling {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 4.GB, task.attempt ) }
    cpus 2
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{*.fasta,*.tsv}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "augur_filter_subsampling.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "augur_filter_subsampling.cfg" }
    input:
      path(sequences)
      path(seq_index)
      path(metadata)
    output:
      path 'subsampled_metadata.tsv', emit: subsampled_metadata
      path 'subsampled_sequences.fasta', emit: subsampled_sequences
      path '{*.sh,*.log}', hidden: true
    script:
      query = optWrap('_augur__filter_subsampling_query', '--query "{}"')    
      group_by = optWrap('augur__filter_subsampling_group_by', '--group-by {}')   
      sequences_per_group = optWrap('augur__filter_subsampling_sequences_per_group', '--sequences-per-group {}')        
      extra = optional('augur__filter_subsampling_extra')
      metadata_id_column = optionalOrDefault('metadata-id-columns', "Id name strain") 
      """
        echo ""
        augur filter \
            --sequences ${sequences} \
            --metadata ${metadata} \
            --sequence-index ${seq_index} \
            --metadata-id-columns ${metadata_id_column} \
           ${query} ${group_by} ${sequences_per_group} ${extra} \
            --output subsampled_sequences.fasta \
            --output-metadata subsampled_metadata.tsv 
      """
}

process augur_index {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 1.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.tsv' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "augur_index.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "augur_index.cfg" }
    input:
      path(sequences)
    output:
      path 'seq_index.tsv', emit: seq_index
      path '{*.sh,*.log}', hidden: true
    script:
      extra = optional('augur__index_extra')  
      """    
        augur index \
        --sequences ${sequences}  \
        --output seq_index.tsv ${extra}
      """
}

process augur_mask {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'masked_alignment.fasta' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "mask.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "mask.cfg" }
    input:
      path(alignment)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'masked_alignment.fasta', emit: masked_alignment
    script:
      bed_content = optional('_augur__mask_bed_content')  
      bed_content_par = optWrap('_augur__mask_bed_content', '--mask mask.bed')
      extra = optional('augur__mask_extra')  
      """
        echo -e "${bed_content}" > mask.bed
        augur mask \
          -s ${alignment} \
          ${bed_content_par} \
          --mask-invalid \
          -o masked_alignment.fasta ${extra}
      """
}