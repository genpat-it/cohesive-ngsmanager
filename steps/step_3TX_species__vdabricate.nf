nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata; taskMemory } from '../functions/common.nf'
include { param;getInput } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def abricateDBDir = param('step_3TX_species__vdabricate__db')

def ex = executionMetadata()

def STEP = '3TX_species'
def METHOD = 'vdabricate' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process abricate {
    container 'quay.io/biocontainers/abricate:0.9.8--h1341992_0'
    containerOptions = " -v ${abricateDBDir}:${abricateDBDir}:ro"
    memory { taskMemory( 6.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(scaffolds), val(database)
    output:
      path '*'
      tuple val(riscd), path("${base}_vdabricate_calls.txt"), emit: calls
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*calls.txt'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_vdabricate.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      riscd = getRisCd(md, ex, STEP, METHOD)
      """                
        abricate ${scaffolds} --datadir ${abricateDBDir} -db ${database} &>> ${base}_vdabricate.log >> ${base}_vdabricate_calls.txt        
      """
}

workflow step_3TX_species__vdabricate {
    take: data
    main:
      abricate(data);
    emit:
      calls = abricate.out.calls
}

workflow {
    getInput()
        .combine(['viruses_TREF'])
        .set { scaffoldsAbricatedb }
    step_3TX_species__vdabricate(scaffoldsAbricatedb)
}
