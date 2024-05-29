nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;param;isCompatibleWithSeqType;isIlluminaPaired   } from '../functions/parameters.nf'
include { stepInputs;parseRISCD } from '../functions/common.nf'

def ex = executionMetadata()

process module_qc_fastqc {
    container 'biocontainers/fastqc:v0.11.5_cv4'
    memory { taskMemory( 1.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    maxForks 10
    when:
      isCompatibleWithSeqType(reads, ['illumina_paired','ion'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true 
    afterScript "echo '${stepInputs(riscd_input, md2, [dt: md2.dt], md2.acc, md2.met, [seq_type:seq_type])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{*.zip,*.html}', saveAs: { filename -> filename.replaceFirst("-DT\\d+_", "-${ex.dt}_") }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      md2 = parseRISCD(riscd_input)   
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      seq_type = isIlluminaPaired(reads) ? 'illumina_paired' : 'ion'
      """
      fastqc $reads &> "${base}_fastqc.log" 
      """
}

workflow {
    module_qc_fastqc(getInput())
}

