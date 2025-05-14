nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;param;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4AN_AMR'
def METHOD = 'resfinder'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process resfinder {
    container "ghcr.io/genpat-it/resfinder:4.1.5--858219071c"
    containerOptions = "-u 0:0"
    memory { taskMemory( 4.GB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
      val genus_species
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.txt,*.fsa,*_format_*.json}', saveAs: { f -> "${base}_${f}" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, '']
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      species = genus_species.replace('_', ' ').toLowerCase()
      """
        run_resfinder.py -o ./ -s "${species}" -l 0.6 -t 0.9 --acquired --db_path_res /resfinder/db_resfinder -ifq ${r1} ${r2} 
      """      
}

workflow step_4AN_AMR__resfinder {
    take: 
      reads
      genus_species
    main:
      resfinder(reads, genus_species);
}

workflow {
    step_4AN_AMR__resfinder(getSingleInput(),param('genus_species'))
}
