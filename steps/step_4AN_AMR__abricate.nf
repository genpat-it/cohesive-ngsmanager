nextflow.enable.dsl=2

include {  flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4AN_AMR'
def METHOD = 'abricate'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def dbs = '"argannot" "card" "ecoh" "ecoli_vf" "ncbi" "plasmidfinder" "resfinder" "vfdb"'

/*
  si potrebbe usare 'each' e passare array di dbs come input
  poi collezionare i risultati con collectFile
  MA bisogna vedere come gestire le esecuzioni con più di un campione (e.g. da filesheet)
  le esecuzioni di abricate sono abbstanza veloci al momento da poter essere eseguite sequenzialmente nello stesso container
*/
process abricate {
    container 'quay.io/biocontainers/abricate:0.9.8--h1341992_0'
    containerOptions = "-v /mnt/biowork:/mnt/biowork:ro -v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    memory { taskMemory( 250.MB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
        input:
      tuple val(riscd_input), path(scaffolds200)
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*_import_abricate.csv,*.json}'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_abricate.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*output_*.csv,*.summary,*calls.txt}'    
    script:
      md = parseMetadataFromFileName(scaffolds200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """                
        for db in ${dbs}
        do
          abricate ${scaffolds200} -db \${db} &>> ${base}_abricate.log >> ${base}_abricate_calls.txt        
        done
        abricate --summary ${base}_abricate_calls.txt > ${base}_abricate.summary
        /scripts/process-abricate-result.py ${base}_abricate.summary ${md.ds} ${md.cmp} ${md.dt}
      """
}

workflow step_4AN_AMR__abricate {
    take: data
    main:
      abricate(data);
}

workflow {
    step_4AN_AMR__abricate(getInput())
}
