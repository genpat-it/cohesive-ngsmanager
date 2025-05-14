nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;isSpeciesSupported;taskMemory } from '../functions/common.nf'
include { getSingleInput;param;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '3QC_coverage'
def METHOD = 'truecoverage' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_SPECIES_ALLOWED = [
  "Escherichia_coli",
  "Salmonella_enterica",
  "Listeria_monocytogenes",
  "Streptococcus_agalactiae",
  "Streptococcus_dysgalactiae",
  "Streptococcus_pneumoniae",
  "Streptococcus_pyogenes",
  "Yersinia_enterocolitica",
  "Haemophilus_influenzae"
]

process truecoverage {
    container "flowcraft/true_coverage:3.2-1"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 8.GB, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, 'illumina_paired', task.process) && isSpeciesSupported(genus_species, GENUS_SPECIES_ALLOWED, reads, task.process)
    input:
      tuple val(riscd_input), path(reads)
      val(genus_species)
    output:
      path '*'
      // path '*_report.csv'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*sample_data_general*.json', saveAs: { "${base}_report.json" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*_report*.txt', saveAs: { "${base}_report.txt" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = reads
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_truecoverage"
      species = genus_species.replace('_', ' ')
      """
        trueCoverage_rematch.py -f ${r1} ${r2} --species ${species} -i /NGStools/true_coverage/data --threads ${params.max_cpus} --json
      """
}

workflow step_3QC_coverage__truecoverage {
    take: 
      reads
      genus_species
    main:
      truecoverage(reads, genus_species)
}

workflow {
    step_3QC_coverage__truecoverage(getSingleInput(),param('genus_species'))
}