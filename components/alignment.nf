nextflow.enable.dsl=2

include { param;optional;paramWrap;optionalOrDefault } from '../functions/parameters.nf'
include { taskCpus;taskMemory } from '../functions/common.nf'

process maaft {
    container "quay.io/biocontainers/mafft:7.520--h031d066_3"
    memory { taskMemory( 150.GB, task.attempt ) }    
    cpus params.max_cpus
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "maaft.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "maaft.cfg" }
    input:
      path(multifasta)
    output:
      path 'alignment.fasta', emit: alignment
      path '{*.sh,*.log}', hidden: true
    script:     
      extra = optionalOrDefault('maaft__extra', '--auto')
      """
        mafft \
        --thread ${task.cpus} \
        ${extra} \
        ${multifasta} \
         > alignment.fasta
      """
}

process gubbins {  
    container "${params.local_registry}/bioinfo/gubbins:3.4--5d60b3c5a2"
    memory { taskMemory( 4.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'gubbins'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "gubbins.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "gubbins.cfg" }
    input:
      path(alignment)
    output:
      path '*'
      path 'gubbins/gubbins_out.masked.aln', emit: masked_alignment
      path 'gubbins/gubbins_out.filtered_polymorphic_sites.fasta', emit: filtered_polymorphic_sites
      path 'gubbins', emit: result
      path '{*.sh,*.log}', hidden: true
    script:
      """
      /usr/local/bin/_entrypoint.sh run_gubbins.py --prefix gubbins_out  ${alignment}
      /usr/local/bin/_entrypoint.sh mask_gubbins_aln.py --aln ${alignment} --gff gubbins_out.recombination_predictions.gff --out gubbins_out.masked.aln
      mkdir gubbins && mv gubbins_out.* gubbins/
      """
}

process snippy_core {
    container "${params.local_registry}/bioinfo/snippy:4.5.1--7be4a1c45a"
    memory { taskMemory( 4.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{core*,*.aln}'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "snippy_core.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "snippy_core.log" }
    input:
      path vcf_folders, stageAs: 'data?'
      tuple val(ref_riscd), val(ref_code), path(ref_file)
    output:
      path '*'
      path 'core.aln', emit: core_alignment
      path 'clean.full.aln', emit: full_alignment
      path '{*.sh,*.log}', hidden: true
    script:
      // XXX renaming folders getting sample name from the first vcf file inside
      """          
        #!/bin/bash -euo pipefail
        mkdir inputs && cd inputs && for dir in ${vcf_folders}; do ln -s ../\${dir} `ls ../\${dir}/*.vcf | head -n 1 | sed -E 's/.+DS[[:digit:]]+-DT[[:digit:]]+_([^_]+)_[[:graph:]]+/\\1/'`; done && cd ..
        snippy-core --ref ${ref_file} --prefix core --inprefix snps inputs/*
        snippy-clean_full_aln core.full.aln > clean.full.aln
      """
}

process muscle {
    container "quay.io/biocontainers/muscle:5.3--h4ac6f70_0"
    memory { taskMemory( 150.GB, task.attempt ) }        
    cpus params.max_cpus
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "muscle.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "muscle.cfg" }
    input:
      path(multifasta)
    output:
      path 'alignment.fasta', emit: alignment
      path '{*.sh,*.log}', hidden: true
    script:
      aln_command = paramWrap('muscle__command', '-{}')
      extra = optional('muscle__extra')
      """
      muscle ${aln_command} \
        ${multifasta} \
        -output alignment.fasta \
        -threads ${task.cpus} \
        ${extra}
      """
}

process prepare_multifasta {
    container "ubuntu:20.04"
    memory { taskMemory( 1.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "prepare_multifasta.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "prepare_multifasta.cfg" }
    input:
      path(summary)
      path(sequences)
    output:
      path 'multifasta.fasta', emit: multifasta
      path '{*.sh,*.log}', hidden: true
    script:
      METADATA_ID_FULL = param('metadata_with_full_id') as boolean    
      """
        #!/bin/bash -euo pipefail
        while IFS=\$'\t' read -r key file; do
            if (${METADATA_ID_FULL}); then
                KEY=\$key
            else
                KEY=\$(sed -E 's/DS[[:digit:]]+-DT[[:digit:]]+_([[:digit:]]+\\.[[:alpha:]]+\\.[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+).+/\\1/' <<< "\$file")
            fi
            awk -v id="\$KEY" '
              NR==1 {\$0=">" id} 
              {print}
            ' \$file
        done < ${summary} > multifasta.fasta       
      """      
}

process amas_summary {
  container 'quay.io/biocontainers/amas:1.0--pyh864c0ab_0'
  memory { taskMemory( 4.GB, task.attempt ) }
  cpus { taskCpus(32, task.attempt) }   
  publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'amas_summary.tsv'
  publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*-seq-summary.txt', saveAs: { "amas_seq_summary.tsv" }
  publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "amas_summary.log" }
  publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "amas_summary.cfg" }
  input:
    path(alignment)
  output:
    path '*'
    path '{*.sh,*.log}', hidden: true
  script:         
    """
      AMAS.py summary -f fasta -d dna -i ${alignment} -c ${task.cpus} -s -e -o amas_summary.tsv
    """
}

process constant_sites {
  container 'quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0'
  memory { taskMemory( 2.GB, task.attempt ) }
  publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'constant_sites', saveAs: { "snpsites_constant_sites.txt" }
  publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "snpsites.log" }
  publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "snpsites.cfg" }
  input:
    path(alignment)
  output:
    path 'constant_sites.txt', emit: constant_sites
    path '{*.sh,*.log}', hidden: true
  script:         
    """
      snp-sites -C ${alignment} > constant_sites.txt
    """
}

process seqkit_exclude {
  container 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'
   memory { taskMemory( 4.GB, task.attempt ) }
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: 'filtered_alignment.fasta', saveAs: { "seqkit_filtered_alignment.fasta" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "seqkit.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "seqkit.cfg" }
    input:
      path(alignment)
      val(sample_to_exclude)
    output:
      path 'filtered_alignment.fasta', emit: filtered_alignment
      path '{*.sh,*.log}', hidden: true
    script:         
      """
        seqkit grep -v -p '"${sample_to_exclude}"' ${alignment} > filtered_alignment.fasta
      """
}

process seqkit_stats {
  container 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'
  memory { taskMemory( 4.GB, task.attempt ) }
  publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.tsv'
  publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "seqkit_stats.log" }
  publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "seqkit_stats.cfg" }
  input:
    path(alignment)
  output:
    path '*'
    path '{*.sh,*.log}', hidden: true
  script:         
    """
      seqkit stats -a -o seqkit_stats.tsv ${alignment}
    """
}