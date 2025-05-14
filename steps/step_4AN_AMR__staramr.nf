nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;isSpeciesSupported } from '../functions/common.nf'
include { getSingleInput;param } from '../functions/parameters.nf'
include { stepInputs;flattenPath } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4AN_AMR'
def METHOD = 'staramr' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_ALLOWED = [
  'campylobacter'
]

def POINFINDER_ORGANISM = [
  campylobacter : "campylobacter"
]

def getPointfinderParam(gsp, map) {
 try {  
    def genus_species = gsp ? gsp.toLowerCase() : ''
    def (genus, species) = genus_species.contains("_") ? genus_species.split('_') : [ genus_species, null ]

    if (map.containsKey(genus)) {
      return map.get(genus)
    } 
    exit 1, "could not find pointfinder param for: ${gsp}"
  } catch(Throwable t) {
      exit 1, "unexpected exception: ${t.asString()}"
  } 
}

process staramr {
    container "ghcr.io/genpat-it/staramr:0.9.1--8fe6b5a239"
    containerOptions = "--user root"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 2.GB, task.attempt ) }
    when:
      isSpeciesSupported(genus_species, GENUS_ALLOWED, assembly, task.process)
    input:
      tuple val(riscd_input), path(assembly)
      val(genus_species)
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'result/{*.tsv,*.xlsx}', saveAs: { filename -> "${base}_${flattenPath(filename)}" }  
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'result/hits/*', saveAs: { filename -> "hits/${base}_hits_${flattenPath(filename) -~ /_DS.+/}.fasta" }  
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: 'result/*.txt', saveAs: { filename -> "${base}_${flattenPath(filename)}" }  
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_staramr"
      pointfinder_organism = getPointfinderParam(genus_species, POINFINDER_ORGANISM)
      """
        sed 's/^>.*/>/g' ${assembly} | awk '{for(i=1;i<=NF;i++){if(\$i~/^>/){\$i=">contig"++count}}} 1' > ${md.cmp}.fasta
        staramr search --pointfinder-organism ${pointfinder_organism} -o result ${md.cmp}.fasta
      """
}

workflow step_4AN_AMR__staramr {
    take: 
      assembly
      genus_species
    main:
      staramr(assembly, genus_species)
}

workflow {
    step_4AN_AMR__staramr(getSingleInput(), param('genus_species'))
}