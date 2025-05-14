nextflow.enable.dsl=2

include { stepInputs;isImportedRiscd;parseRISCD;parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;param  } from '../functions/parameters.nf'

def ex = executionMetadata()

process quast {
    container 'quay.io/biocontainers/quast:4.4--boost1.61_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }   
    maxForks 10
    input:
      tuple val(riscd_input), path(assembly)
    output:
      path '{*_quast.tsv,*.json}', emit: tsv
      path '{*.sh,*.log}', hidden: true 
    afterScript "echo '${stepInputs(riscd_input, md_input, [dt: md_input.dt], md_input.acc, md_input.met, null)}' > ${base}_input.json${!imported ? '.ignore' : ''}"
    publishDir mode: 'rellink', "${params.outdir}${result_path}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}${result_path}/qc/result", pattern: '*.tsv'   
    publishDir mode: 'rellink', "${params.outdir}${result_path}/qc/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}${result_path}/qc/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
        md = parseMetadataFromFileName(assembly.getName())
        md_input = parseRISCD(riscd_input)    
        imported = isImportedRiscd(riscd_input)   
        result_path = imported ? "/${md.anno}/${md.cmp}/${md_input.acc}/${md_input.ds}-${md_input.dt}_${md_input.met}" : ''
        base = "${md.ds}-${ex.dt}_${md.cmp}_quast"
        """
        quast -m 200 --fast -o quast ${assembly}
        cut -f1,14,15,16,17,18,19,20,21 quast/transposed_report.tsv > ${base}.tsv 
        """
}

process summary {
    container 'ubuntu:20.04'
    memory { taskMemory( 200.MB, task.attempt ) }
    when:
      (quast_results instanceof java.util.Collection) && quast_results.size() > 1      
    input:
      path(quast_results)
    output:
      path '*.tsv'
      path '{*.sh,*.log}', hidden: true 
    publishDir mode: 'rellink', "${params.outdir}/result", pattern: '*.tsv'   
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
        base = "summary"
        """
        cat ${quast_results} | sort -u > ${base}.tsv
        """
}

workflow module_qc_quast {
    take: 
      input
    main:
      quast_result = quast(input).tsv
      if (params.module_qc_quast__summary) {
        summary(quast_result.collect())
      }
}

workflow {
    tsv = quast(getInput()).tsv.collect()
    if (params.module_qc_quast__summary) {
      summary(tsv)
    }
}

