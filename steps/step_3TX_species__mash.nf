nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;isCompatibleWithSeqType } from '../functions/parameters.nf'

def DB_PATH="/mnt/biowork/databases/PROGRAMS/mash"

def ex = executionMetadata()

def STEP = '3TX_species'
def METHOD = 'mash' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process mash {
    container "${LOCAL_REGISTRY}/bioinfo/mash:2.3--debdd7eb23"
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
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, '']
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_mash"
      """
        mash screen -w -p ${task.cpus} /mash/refseq/refseq.genomes.k21s1000.msh ${r1} ${r2} > ${base}_screen.tsv
        set +e  # allow continuing even if nonzero exit code
        assess_mash_screen.py --gram /mash/BACTpipe/resources/gram_stain.txt ${base}_screen.tsv > ${base}_screening_results.tsv
        exit_code=$?
        if [ "$exit_code" -eq 0 ]; then
            echo "Screening succeeded."
        elif [ "$exit_code" -eq 3 ]; then
            echo "Contamination detected, proceeding."
        else
            echo "Screening failed with exit code $exit_code."
            exit "$exit_code"
        fi     
      """
}

workflow step_3TX_species__mash {
    take: 
      reads
    main:
      mash(reads)
}

workflow {
    step_3TX_species__mash(getInput())
}