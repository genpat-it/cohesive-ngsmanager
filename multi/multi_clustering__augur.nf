nextflow.enable.dsl=2

include { taskMemory;getEmpty } from '../functions/common.nf'
include { optional;optWrap;getInput;getReference;param;optionalOrDefault;optionalBoolean } from '../functions/parameters.nf'
include { sequence_merge;metadata_merge;geodata_merge } from '../components/misc.nf'
include { prepare_multifasta;amas_summary;seqkit_stats } from '../components/alignment.nf'
include { augur_filter as filter;augur_filter_subsampling as subsampling;augur_index as index;augur_mask as mask } from '../components/augur.nf'
include { reportree_alignment;prepare_reportree_metadata } from '../components/clustering.nf'

process align {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 15.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'alignment.fasta' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "align.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "align.cfg" }
    input:
      path(sequences)
      tuple val(ref_sample), val(ref_resultcode), val(ref_code), path(ref_path)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'alignment.fasta', emit: alignment
    script:
      extra = optional('augur__align_extra')  
      remove_reference = optionalBoolean('augur__align_remove_reference') ?  '--remove-reference' : ''
      """     
        if grep -q "^>${ref_sample}\$" ${sequences}; then
            ref_option="--reference-name ${ref_sample}"
        else
            ref_option="--reference-sequence ${ref_path}"
        fi
        echo "ref_option=\${ref_option}"
        augur align \
            --nthreads auto \
            --sequences ${sequences} \
            \$ref_option \
            --output alignment.fasta \
            --fill-gaps ${remove_reference} ${extra}
      """
}

process tree {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'tree_raw.nwk' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "tree.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "tree.cfg" }
    input:
      path(alignment)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'tree_raw.nwk', emit: tree_raw
    script:
      method = param('augur__tree_method')  
      extra = optional('augur__tree_extra')  
      """
        augur tree \
          --alignment ${alignment} \
          --output tree_raw.nwk \
          --method ${method} \
          --nthreads auto ${extra}
      """
}

process refine {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 2.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'tree.nwk' 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'branch_lengths.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "refine.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "refine.cfg" }
    input:
      path(tree)
      path(alignment)
      tuple val(ref_sample), val(ref_resultcode), val(ref_code), path(ref_path)
      path(metadata)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'tree.nwk', emit: tree
      path 'branch_lengths.json', emit: branch_lengths
    script:
      coalescent = param('augur__refine_coalescent')  
      reroot = optionalOrDefault('augur__refine_reroot', 'best')
      remove_reference = optionalBoolean('augur__align_remove_reference') 
      extra = optional('augur__refine_extra')  
      """
        if [ "${reroot}" == "reference" ] && [ "${remove_reference}" == "false" ]; then
          if [ "`cut -f1 ${metadata} | grep -c '^${ref_sample}'`" == "1" ]; then
              ROOTING='--root ${ref_sample}'
          else
              ROOTING='--root ${ref_code}'
          fi
        elif [ "${reroot}" == "reference" ] && [ "${remove_reference}" == "true" ]; then
          ROOTING='--root best'
        else
          ROOTING='--root ${reroot}'
        fi
        echo "ROOTING=\${ROOTING}"
         augur refine \
            --tree ${tree} \
            --alignment ${alignment} \
            --metadata ${metadata} \
            --output-tree tree.nwk \
            --output-node-data branch_lengths.json \
            --timetree \
            --coalescent ${coalescent} \
            --date-confidence \$ROOTING ${extra}
      """
}

process ancestral {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 2.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'nt_muts.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "ancestral.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "ancestral.cfg" }
    input:
      path(tree)
      path(alignment)
      tuple val(ref_sample), val(ref_resultcode), val(ref_code), path(ref_path)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'nt_muts.json', emit: nt_muts
    script:
      inference = param('augur__ancestral_inference')  
      extra = optional('augur__ancestral_extra')      
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
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 2.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'aa_muts.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "translate.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "translate.cfg" }
    input:
      path(tree)
      path(nt_muts)
      tuple val(ref_sample), val(ref_resultcode), val(ref_code), path(ref_path)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'aa_muts.json', emit: aa_muts
    script:
      extra = optional('augur__translate_extra')      
      """
        augur translate \
            --tree ${tree} \
            --ancestral-sequences ${nt_muts} \
            --reference-sequence ${ref_path} \
            --output aa_muts.json ${extra}
      """
}

process traits {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 2.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'traits.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "traits.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "traits.cfg" }
    input:
      path(tree)
      path(metadata)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'traits.json', emit: traits
    script:
      columns = param('augur__traits_columns')      
      extra = optional('augur__traits_extra')      
      """
        augur traits \
            --tree ${tree} \
            --metadata ${metadata} \
            --output-node-data traits.json \
            --columns ${columns} \
            --confidence ${extra}
      """
}


process clades {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 2.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'clades.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "clades.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "clades.cfg" }
    input:
      path(tree)
      path(nt_muts)
      path(aa_muts)
      path(clades)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'clades.json', emit: clades
    script:
      clades_content = optional('_augur__clades_content')
      extra = optional('augur__clades_extra')      
      if (clades.exists() && !clades.empty()) {
        """
          augur clades \
              --tree ${tree} \
              --mutations ${nt_muts} ${aa_muts} \
              --clades ${clades} \
              --output-node-data clades.json ${extra}
        """
      } else if (clades_content) {
        """
          echo -e "${clades_content}" > _clades_definition_.tsv
          augur clades \
              --tree ${tree} \
              --mutations ${nt_muts} ${aa_muts} \
              --clades _clades_definition_.tsv \
              --output-node-data clades.json ${extra}
        """
      } else {
        """
          echo "{}" > clades.json
        """
      }

}

process prepare_metadata {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 1.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'augur_metadata.tsv' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "prepare_metadata.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "prepare_metadata.cfg" }
    input:
      path(metadata)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'augur_metadata.tsv', emit: metadata
    script:
      ID_COLUMN = param('metadata_id_column')
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
        ' ${metadata} | sed 's/${ID_COLUMN}/name/i' | sed -E 's/${SUMMARY_DATE_ALIASES}/date/i' > augur_metadata.tsv	
      """
}

process filter_by_partition {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    containerOptions = "-v ${workflow.projectDir}/scripts/module_augur:/scripts:ro"
    memory { taskMemory( 15.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'selection.txt', saveAs: { "samples_selected_for_tree_building.txt" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "filter_by_partition.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "filter_by_partition.cfg" }
    input:
      path(summary)
      path(alignment)
      path(metadata)
    output:
      path '{*.sh,*.log}', hidden: true
      path 'sequences_filter_by_partition.fasta', emit: alignment
      path 'metadata_filter_by_partition.tsv', emit: metadata
    script:
      partition = optional('filter_before_tree_by_partition')  
      metadata_id_column = optionalOrDefault('metadata-id-columns', "Id name strain") 
      """     
        /scripts/filter_by_partition.py --partition_summary ${summary} --partition "MST-${partition}x1.0" --output_samples selection.txt

        augur filter --exclude-all --include selection.txt --metadata ${metadata} --sequences ${alignment} --output-metadata metadata_filter_by_partition.tsv --output-sequences sequences_filter_by_partition.fasta --metadata-id-columns ${metadata_id_column} 
      """
}

process filter_by_closest {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    containerOptions = "-v ${workflow.projectDir}/scripts/module_augur:/scripts:ro"
    memory { taskMemory( 15.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'selection.txt', saveAs: { "samples_selected_for_tree_building.txt" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "filter_by_closest.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "filter_by_closest.cfg" }
    input:
      path(matrix)
      path(alignment)
      path(metadata)
    output:
      path '{*.sh,*.log}', hidden: true
      path 'sequences_filter_by_closest.fasta', emit: alignment
      path 'metadata_filter_by_closest.tsv', emit: metadata
    script:
      closest_num = optional('filter_before_tree_by_closest')  
      soi = optional('reportree__sample_of_interest').replaceAll(/[\s\n\t\r"'$\{]/, "")	
      metadata_id_column = optionalOrDefault('metadata-id-columns', "Id name strain") 
      """
        /scripts/filter_by_closest.py --include_input --matrix ${matrix} --samples "${soi}" --closest_num ${closest_num} --output selection.txt

        augur filter --exclude-all --include selection.txt --metadata ${metadata} --sequences ${alignment} --output-metadata metadata_filter_by_closest.tsv --output-sequences sequences_filter_by_closest.fasta --metadata-id-columns ${metadata_id_column} 
      """
}

process export {
    container "quay.io/biocontainers/augur:29.0.0--pyhdfd78af_0"
    memory { taskMemory( 2.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'auspice.json' 
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "export.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "export.cfg" }
    input:
      path(tree)
      path(metadata)
      path(branch_lengths)
      path(traits)
      path(nt_muts)
      path(aa_muts)
      path(clades)
      path(lat_longs)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path 'auspice.json', emit: auspice
    script:
      extra = optional('augur__export_extra')  
      //TODO --auspice-config {input.auspice_config} \
      GEO_RESOLUTION_COLUMNS = param('metadata_geo_column')      
      """
        METADATA_LIST=\$(head -n 1 ${metadata} | tr \$'\t' ' ')
        augur export v2 \
            --tree ${tree} \
            --metadata ${metadata} \
            --node-data ${branch_lengths} ${traits} ${clades} ${nt_muts} ${aa_muts} \
            --lat-longs ${lat_longs} \
            --color-by-metadata \${METADATA_LIST} \
            --geo-resolutions ${GEO_RESOLUTION_COLUMNS} \
            --output auspice.json ${extra}
      """
}

workflow multi_clustering__augur {
    take: 
      input
      reference
      raw_metadata
      _geodata
      clades_definition
    main:  
      summary = input.collectFile { [ "summary.tsv", it[0] + '\t' + (it[1][0]?.getName() ?: '') + '\n' ] }
      files= input.collect { it[1] }
      
      _multifasta = prepare_multifasta(summary, files).multifasta
      ext_seq = optional('ext_sequences')
      ext_metadata = optional('ext_metadata')
      ext_geodata = optional('ext_geodata')

      // get first available metadata - if possibile, use extra columns added by reportree
      if (ext_seq && ext_metadata && ext_geodata) {
        multifasta = sequence_merge(_multifasta, ext_seq).merged_sequences
        _metadata = metadata_merge(raw_metadata, ext_metadata).merged_metadata
        geodata = geodata_merge(_geodata, ext_geodata).merged_geodata
      } else {
        multifasta = _multifasta
        _metadata = raw_metadata
        geodata = _geodata
      }
      metadata = prepare_metadata(_metadata).metadata
      seq_index = index(multifasta).seq_index
      filtered_metadata = filter(multifasta, seq_index, metadata).filtered_metadata
      seqkit_stats(filter.out.filtered_sequences)
      raw_alignment = align(filter.out.filtered_sequences, reference).alignment
      alignment = mask(raw_alignment).masked_alignment
      if (optionalBoolean('ma_summary')) {
        amas_summary(alignment)
      }
      if (optional('reportree__thr')) {
        reportree_alignment(alignment, filtered_metadata,  'grapetree', '/tmp/fake')      
        if (optional('reportree__sample_of_interest')) {
          if (optional('filter_before_tree_by_partition')) {
            filter_by_partition(reportree_alignment.out.soi_partition_summary, alignment, filtered_metadata)
            metadata_pre_subsampling = filter_by_partition.out.metadata
            alignment_pre_subsampling = filter_by_partition.out.alignment
            if (optional('filter_before_tree_by_closest')) {
              log.warn "'filter_before_tree_by_partition' is specified, so 'filter_before_tree_by_closest' will be ignored (these parameters cannot be used together)"
            }
          } else if (optional('filter_before_tree_by_closest')) {
            filter_by_closest(reportree_alignment.out.matrix, alignment, filtered_metadata)
            metadata_pre_subsampling = filter_by_closest.out.metadata
            alignment_pre_subsampling = filter_by_closest.out.alignment
          } else {
            metadata_pre_subsampling = reportree_alignment.out.metadata_w_partitions
            alignment_pre_subsampling = alignment
          }
        } else {
          if (optional('filter_before_tree_by_closest') || optional('filter_before_tree_by_partition') ) {
            log.warn "'reportree__sample_of_interest' not specified, so 'filter_before_tree_by_closest' and 'filter_before_tree_by_partition' will be ignored"
          }
          metadata_pre_subsampling = reportree_alignment.out.metadata_w_partitions
          alignment_pre_subsampling = alignment
        }           
      } else {
        if (optional('filter_before_tree_by_closest') || optional('filter_before_tree_by_partition') ) {
          log.warn "'reportree__thr' not specified, so 'filter_before_tree_by_closest' and 'filter_before_tree_by_partition' will be ignored"
        }
        metadata_pre_subsampling = filtered_metadata
        alignment_pre_subsampling = alignment
      }

      metadata_post_subsampling = subsampling(alignment_pre_subsampling, seq_index, metadata_pre_subsampling).subsampled_metadata
      alignment_post_subsampling = subsampling.out.subsampled_sequences

      tree_raw = tree(alignment_post_subsampling).tree_raw
      tree = refine(tree_raw, alignment_post_subsampling, reference, metadata_post_subsampling).tree
      nt_muts = ancestral(tree, alignment_post_subsampling, reference).nt_muts
      aa_muts = translate(tree, nt_muts, reference).aa_muts
      traits = traits(tree, metadata_post_subsampling).traits
      clades = clades(tree, nt_muts, aa_muts, clades_definition ).clades
      export(tree, metadata_post_subsampling, refine.out.branch_lengths, traits, nt_muts, aa_muts, clades, geodata)
  }

workflow {
    reference = getReference('gb')    
    multi_clustering__augur(getInput(), reference, param('metadata'), param('geodata'), optionalOrDefault('clades', getEmpty()))
}
