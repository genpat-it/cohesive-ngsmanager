nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;isCompatibleWithSeqType;param } from '../functions/parameters.nf'

//def DB_PATH="/mnt/biowork/databases/PROGRAMS/mash"
//def MSH_DATASET= "/mnt/biowork/home/p.castelli/nf_pipeline_test-20240715/cansnp-mash"

def ex = executionMetadata()

def STEP = '4TY_distance'
def METHOD = 'mash' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def DB_PATH = param("step_4TY_distance__mash__db_path")
def DB_NAME = param("step_4TY_distance__mash__db_name")

process mash_sketch {
    container "${LOCAL_REGISTRY}/bioinfo/mash:2.3--debdd7eb23"
    containerOptions = "-v ${DB_PATH}:/dataset_msh:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 10.GB, task.attempt ) }
    cpus 8
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.tsv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.msh'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, '']
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_mash"
        """
          zcat ${r1} ${r2} | mash sketch -r -m 2 -o ${base} -
        """ 
//          mash dist /dataset_msh/${DB_NAME} ${base}.msh > ${md.ds}-${ex.dt}_${md.cmp}_dist.tsv
}

workflow step_4TY_distance__mash {
    take: 
      reads
    main:
      mash_sketch(reads)
}

workflow {
    step_4TY_distance__mash(getInput())
}