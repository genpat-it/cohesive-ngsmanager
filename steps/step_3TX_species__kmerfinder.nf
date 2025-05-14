nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata; csv2map; extractKey;taskMemory } from '../functions/common.nf'
include { param;getInput } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

KMERFINDER_SPECIES_DIR = param('step_3TX_species__kmerfinder__db')
KMERFINDER_REFERENCE_DIR = "${KMERFINDER_SPECIES_DIR}/Bacteria/Fasta/"

def ex = executionMetadata()

def STEP = '3TX_species'
def METHOD = 'kmerfinder' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def getBacterialReferencePath(checkFile) {
    try {
      bacterialReferenceId = csv2map(checkFile, "\\t").assembly_accBacteria
      return [ '-', bacterialReferenceId, file("${KMERFINDER_REFERENCE_DIR}/${bacterialReferenceId}*.fa") ] // [riscd (empty), refCode, refPath]
    } catch(Throwable t) {
        exit 1, "could not get bacterialReferenceId from '${checkFile}', exception: ${t.asString()}"
    }
}

def getCalculatedSpecies(checkFile) {
    try {
      return csv2map(checkFile, "\\t").speciesAssigned
    } catch(Throwable t) {
        exit 1, "could not get speciesAssigned from '${checkFile}', exception: ${t.asString()}"
    }
}

process kmerfinder {
    container "python:2.7.8"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro -v ${KMERFINDER_SPECIES_DIR}:${KMERFINDER_SPECIES_DIR}:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 3.GB, task.attempt ) }
    input:
      tuple val(riscd_input), path(scaffolds200)
    output:
      tuple val(riscd), path("${base}_kmerfinder.check"), emit: check
      path "**"
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*/*.tsv', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.json,*/*.check,*/*.log}', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_kmerfinder.cfg" }
    script:
      md = parseMetadataFromFileName(scaffolds200.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      inputFile = scaffolds200
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
        /scripts/FindTemplate.py -i ${inputFile} -t ${KMERFINDER_SPECIES_DIR}/Bacteria_DB -o ${base}_kmerfinder_bacterial.tsv -x ATGAC -w  > ${base}_kmerfinder.log
        /scripts/FindTemplate.py -i ${inputFile} -t ${KMERFINDER_SPECIES_DIR}/Viral_DB -o ${base}_kmerfinder_viral.tsv -x ATGAC -w  >> ${base}_kmerfinder.log
        /scripts/ParseSpeciesFile_newDB.py ${base}_kmerfinder_bacterial.tsv ${KMERFINDER_SPECIES_DIR}
      """
}

workflow step_3TX_species__kmerfinder {
    take: data
    main:
      kmerfinder(data);
    emit:
      assigned_species = kmerfinder.out.check.map { [ it[0], getCalculatedSpecies(it[1]), getBacterialReferencePath(it[1]) ] }
}

workflow {
    step_3TX_species__kmerfinder(getInput())
}