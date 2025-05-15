nextflow.enable.dsl=2

include { param;optional } from '../functions/parameters.nf'
include { prepare_reportree_metadata;find_closest;fasttree;iqtree;augur;augur as augur_hc;reportree_alignment as reportree_gt; reportree_alignment as reportree_hc} from './clustering.nf'

workflow multi_alignment_analysis {
  take: 
    ma
    raw_metadata
    geodata
    nomenclature
    polymorphic_sites_only         
    fconst     
  main:

    metadata = prepare_reportree_metadata(raw_metadata).metadata

    if (optional('reportree__thr')) {
      reportree_gt(ma, metadata,  'grapetree', nomenclature)
      zoom_gt = reportree_gt.out.zooms
      matrix_gt = reportree_gt.out.matrix
      metadata_w_partitions_gt = reportree_gt.out.metadata_w_partitions
    } else {
      zoom_gt = Channel.empty()
      matrix_gt = Channel.empty()
      metadata_w_partitions_gt = Channel.empty()
    }
    if (optional('reportree__HC_method') || optional('reportree__HC_threshold') || optional('reportree__pct_HC_threshold')) {
      reportree_hc(ma, metadata, 'HC', nomenclature)
      reportree_hc.out.nwk
          .flatten()
          .count()
          .filter { it > 1 }
          .subscribe { log.warn "Reportree-HC generated multiple dendrograms; only one will be displayed in Auspice" }
      metadata_w_partitions_hc = reportree_hc.out.metadata_w_partitions      
      zoom_hc = reportree_hc.out.zooms
      matrix_hc = reportree_hc.out.matrix
      augur_hc(reportree_hc.out.nwk.flatten().first(), metadata_w_partitions_hc, geodata, 'reportree')
    } else {
      zoom_hc = Channel.empty()
      matrix_hc = Channel.empty()
      metadata_w_partitions_hc = Channel.empty()
    }
    if (optional('closest_samples_report_threshold') && optional('reportree__sample_of_interest')) {
      matrix = matrix_gt.concat(matrix_hc).first()
      find_closest(matrix)
    }
    zoom_gt.concat(zoom_hc)
      .splitText()
        .ifEmpty { "" }
          .collectFile(name: "${params.outdir}/result/zooms.txt" )

    method = optional('mltree_method')
    // get first available metadata - if possibile, use extra columns added by reportree
    reportree_metadata = metadata_w_partitions_gt.concat(metadata_w_partitions_hc).concat(metadata).first()
    if (method == 'fasttree') {
      nwk = fasttree(ma).nwk
      augur(nwk, reportree_metadata, geodata, 'fasttree')
    } else if (method == 'iqtree') {
      nwk = iqtree(ma, polymorphic_sites_only, fconst).nwk
      augur(nwk, reportree_metadata, geodata, 'iqtree')
    }          
}