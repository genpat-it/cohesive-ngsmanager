nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;param;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_plasmid'
def METHOD = 'mobsuite'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process mobsuite {
    container "${LOCAL_REGISTRY}/bioinfo/mobsuite:3.1.4--c2c728b1a2"
    memory { taskMemory( 4.GB, task.attempt ) }
    cpus { [32, params.max_cpus as int].min() }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(assembly)
    output:
      path '**'
      tuple val(riscd), path("output/*.fasta"), emit: plasmids, optional: true
      path '{*.sh,*.log}', hidden: true 
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'output/*', saveAs: { f -> "${base}_${flattenPath(f)}" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
      /usr/local/bin/_entrypoint.sh mob_recon -s ${base} -o output -i ${assembly} -n ${task.cpus} --force --unicycler_contigs --run_overhang
      """      
}

workflow step_4TY_plasmid__mobsuite {
    take: 
      reads
    main:
      plasmids = mobsuite(reads).plasmids
    emit:
      plasmids = plasmids
}

workflow {
    step_4TY_plasmid__mobsuite(getSingleInput())
}
