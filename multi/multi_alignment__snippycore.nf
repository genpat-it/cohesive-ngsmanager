nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;getReferenceUnkeyed;getInputFolders;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent  } from '../functions/parameters.nf'

def ex = executionMetadata()

def ENTRYPOINT = "multi_alignment__snippycore"

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
      path vcf_files, stageAs: 'data?'
      tuple val(ref_riscd), val(ref_code), path(ref_file)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{core*}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "snippy_core.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '*.log', saveAs: { "snippy_core.log" }
    script:
      // XXX renaming folders getting sample name from the first vcf file inside
      """          
        #!/bin/bash -euo pipefail
        mkdir inputs && cd inputs && for dir in ${vcf_files}; do ln -s ../\${dir} `ls ../\${dir}/*.vcf | head -n 1 | sed -E 's/.+DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_[[:graph:]]+/\\1/'`; done && cd ..
        snippy-core --ref ${ref_file} --prefix core --inprefix snps inputs/*
      """
}

workflow multi_alignment__snippycore {
    take: 
        reads      
        reference
    main:
        reads.combine(reference)
                .multiMap { 
                    reads: it[0..1] // riscd, R[]
                    reference:  it[2..4] // riscd, code, path
                }.set { input }
        folders = snippy(input.reads,input.reference).results
        snippy_core(folders.collect(), reference)
}

workflow multi_alignment__snippycore_vcf {
    take: 
        input
        reference
    main:
        input
            .map { it[1] }
            .collect()
            .set { vcfs }       
        snippy_core(vcfs, reference)
}

workflow {
    // 1PP_* => snippy + snippycore
    multi_alignment__snippycore(getInput().filter( ~/^.*\/1PP_.*/ ), getReferenceUnkeyed('gb'))
    // 2AS_* => snippycore only
    multi_alignment__snippycore_vcf(getInputFolders().filter( ~/^.*\/2AS_.*/ ), getReferenceUnkeyed('gb'))
}