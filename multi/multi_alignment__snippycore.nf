nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { param;optionalBoolean;optionalOrDefault;optional;getInputOf;getReferenceUnkeyed;getInputFolders;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent  } from '../functions/parameters.nf'
include { gubbins;snippy_core;seqkit_exclude;seqkit_exclude as seqkit_exclude_from_full;constant_sites;amas_summary } from '../components/alignment.nf'
include { multi_alignment_analysis } from '../components/subworkflows.nf'
include { augur_mask as mask } from '../components/augur.nf'

def ex = executionMetadata()

process snippy {
    container "ghcr.io/genpat-it/snippy:4.5.1--7be4a1c45a"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 8.GB, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base_ref}.log" }    
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(reference_path)
    output:
      path 'snippy', emit: results
      path '{*.sh,*.log}', hidden: true
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_snippy_${reference}"
      is_fasta = r1.getName() ==~ /.+\.fa(sta)?$/
      extra_files_cleanup = optionalBoolean('snippy__keep_intermediate_files') ? "" : "-o -name '*.bam'"
      if (is_fasta) {
        """
        trap "find -name '*.sam' ${extra_files_cleanup} -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --ctgs ${r1} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        """     
      } else if (isIlluminaPaired(reads)) {
        """
        trap "find -name '*.sam' ${extra_files_cleanup} -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --R1 ${r1} --R2 ${r2} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        """
      } else if (isIonTorrent(reads)) {
        """
        trap "find -name '*.sam' ${extra_files_cleanup} -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --se ${r1} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        """      
      }      
}

workflow snippycore {
    take: 
        reads      
        pre_calculated_result_folders      
        reference
        metadata
        geodata
        reportree_nomenclature        
    main:
        reads.combine(reference)
                .multiMap { 
                    reads: it[0..1] // riscd, R[]
                    reference:  it[2..4] // riscd, code, path
                }.set { input }
        folders = snippy(input.reads,input.reference).results
        
        vcf_folders = folders.mix(pre_calculated_result_folders.map { it[1] }).collect()

        snippy_core(vcf_folders, reference)
        
        if (param('snippycore__alignment_type') == 'core') {
          if (optionalBoolean('gubbins__remove_recombinations')) {
            log.warn "Proceeding with the full alignment for Gubbins; note that the output filtered polymorphic sites do not represent the strict core genome alignment."
            unfiltered_core_ma = gubbins(snippy_core.out.full_alignment).filtered_polymorphic_sites
            unfiltered_full_ma = gubbins.out.masked_alignment
          } else {
            unfiltered_core_ma = snippy_core.out.core_alignment
            unfiltered_full_ma = snippy_core.out.full_alignment
          }
        
          if (optionalBoolean('exclude_reference')) {
            ma = seqkit_exclude(unfiltered_core_ma, 'Reference').filtered_alignment
            seqkit_exclude_from_full(unfiltered_full_ma, 'Reference')
            // calculate over full alignment, but excluding reference
            fconst = constant_sites(seqkit_exclude_from_full.out.filtered_alignment).constant_sites      
          } else {
            ma = unfiltered_core_ma
            // calculate over full alignment, with no exclusions
            fconst = constant_sites(unfiltered_full_ma).constant_sites      
          }
        } else {
          if (optionalBoolean('gubbins__remove_recombinations')) {
            unfiltered_ma = gubbins(snippy_core.out.full_alignment).masked_alignment
          } else {
            unfiltered_ma = snippy_core.out.full_alignment
          } 
          if (optionalBoolean('exclude_reference')) {
            ma = seqkit_exclude(unfiltered_ma, 'Reference').filtered_alignment
          } else {
            ma = unfiltered_ma
          }
          // calculate over (eventually filtered) alignment alignment
          fconst = constant_sites(ma).constant_sites
        }

        if (optional('augur__mask_bed_content')) {
          final_alignment = mask(ma).masked_alignment
        } else {
          final_alignment = ma
        }
        if (optionalBoolean('ma_summary')) {
          amas_summary(final_alignment)
        }

        multi_alignment_analysis(ma, metadata, geodata, reportree_nomenclature, param('snippycore__alignment_type') == 'core', fconst)
}

workflow {
    // use getInputOf instead of getInput to avoid warnings
    reads = getInputOf('{*.fasta,*.vcf,*.f*q*.gz}').filter( ~/.*\/(1PP_|2AS_denovo\/|2AS_import\/).*/ )
    pre_calculated_result_folders = getInputFolders().filter( ~/.*\/2AS_mapping.*/ )

    // pre_calculated_result_folders have to be calculated using the same reference and format (gb or fasta)
    snippycore(reads, pre_calculated_result_folders, getReferenceUnkeyed('any'), param('metadata'), param('geodata'), optionalOrDefault('reportree__nomenclature', '/tmp/fake'))
}