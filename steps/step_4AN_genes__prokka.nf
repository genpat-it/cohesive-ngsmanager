nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory;extractKey } from '../functions/common.nf'
include { getInput;getKingdom;getReferenceOptional;checkEnum } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4AN_genes'
def METHOD = 'prokka'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

enum KINGDOM {
	Viruses, Bacteria, Archaea, Mitochondria
}

process prokka {
    container 'quay.io/biocontainers/prokka:1.14.5--pl526_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    when:
      checkEnum(KINGDOM, kingdom)
    input:
      tuple val(riscd_input), path(scaffolds200), val(kingdom), val(riscd_ref), val(reference), path(gb)
    output:
      path '**'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_ref], md, ex, STEP, METHOD, [reference:reference, kingdom: kingdom])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'result/*', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: "{*.log,*.json}", saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_prokka.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"   
      if (gb.empty()) {
        baseMethod =  "${base}_prokka"
        gbParam = ""
      } else {
        baseMethod =  "${base}_prokka_${reference}"
        gbParam = "--proteins ${gb}"
      }
      extraOptions = "--centre X --compliant"
      """
        prokka --kingdom ${kingdom} ${extraOptions} --outdir "result" ${gbParam} --prefix ${baseMethod}_result ${scaffolds200} 2> ${baseMethod}.log
      """
}

workflow step_4AN_genes__prokka {
    take: 
      data
    main:
      prokka(data);
}

workflow {
  getInput()
      .cross(getKingdom()) { extractKey(it) }
      .cross(getReferenceOptional('gb')) { extractKey(it) }
      .map { it.flatten() }  // [ riscd assembly ds kingdom ds riscd_ref refid refpath]
      .map { 
          [ it[0], it[1], it[3], it[5], it[6], it[7] ]
      }
      .set { input }
  step_4AN_genes__prokka(input)
}
