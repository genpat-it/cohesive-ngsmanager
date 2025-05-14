nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;param;isCompatibleWithSeqType;isIlluminaPaired   } from '../functions/parameters.nf'
include { stepInputs;parseRISCD } from '../functions/common.nf'

def ex = executionMetadata()

process module_qc_nanoplot {
    container 'quay.io/biocontainers/nanoplot:1.44.1--pyhdfd78af_0'
    memory { taskMemory( 2.GB, task.attempt ) }
    cpus 8
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    maxForks 10
    when:
      isCompatibleWithSeqType(reads, ['nanopore'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true 
    afterScript "echo '${stepInputs(riscd_input, md2, [dt: md2.dt], md2.acc, md2.met, [seq_type:'nanopore'])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{*.txt,*.html}', saveAs: { filename -> filename.replaceFirst("-DT\\d+_", "-${ex.dt}_") }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }

    script:
      md = parseMetadataFromFileName(reads.getName())
      md2 = parseRISCD(riscd_input)       
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      NanoPlot -t ${task.cpus} --fastq ${reads} --tsv_stats --no_static -o . -p ${base}_
      """
}

workflow {
    module_qc_nanoplot(getInput())
}

