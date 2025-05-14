nextflow.enable.dsl=2

include { taskMemory } from '../functions/common.nf'
include { param;getInput;optional;optionalOrDefault;optionalBoolean } from '../functions/parameters.nf'
include { maaft;gubbins;prepare_multifasta;constant_sites;amas_summary;seqkit_stats } from '../components/alignment.nf'
include { multi_alignment_analysis } from '../components/subworkflows.nf'
include { sequence_merge;metadata_merge;geodata_merge } from '../components/misc.nf'
include { augur_filter;augur_index } from '../components/augur.nf'

workflow multi_clustering__maaft {
    take: 
        input
        _metadata
        _geodata
        nomenclature
    main:
        summary = input.collectFile { [ "summary.tsv", it[0] + '\t' + it[1][0].getName() + '\n' ] }
        files= input.collect { it[1] }

        _multifasta = prepare_multifasta(summary, files).multifasta
        ext_seq = optional('ext_sequences')
        ext_metadata = optional('ext_metadata')
        ext_geodata = optional('ext_geodata')

        // get first available metadata - if possibile, use extra columns added by reportree
        if (ext_seq && ext_metadata && ext_geodata) {
          multifasta = sequence_merge(_multifasta, ext_seq).merged_sequences
          metadata = metadata_merge(_metadata, ext_metadata).merged_metadata
          geodata = geodata_merge(_geodata, ext_geodata).merged_geodata
        } else {
          multifasta = _multifasta
          metadata = _metadata
          geodata = _geodata
        }

        seq_index = augur_index(multifasta).seq_index
        augur_filter(multifasta, seq_index, metadata)
        seqkit_stats(augur_filter.out.filtered_sequences)
        ma = maaft(augur_filter.out.filtered_sequences).alignment
        if (optionalBoolean('ma_summary')) {
          amas_summary(ma)
        }
        fconst = constant_sites(ma).constant_sites

        multi_alignment_analysis(ma, augur_filter.out.filtered_metadata, geodata, nomenclature, false, fconst)
}

workflow {
  multi_clustering__maaft(getInput(),  param('metadata'), param('geodata'), optionalOrDefault('reportree__nomenclature', '/tmp/fake'));
}