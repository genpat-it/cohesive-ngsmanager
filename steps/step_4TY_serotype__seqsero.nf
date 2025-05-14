nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;isSpeciesSupported;taskCpus } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isIonTorrent;param } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_serotype'
def METHOD = 'seqsero' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_SPECIES_ALLOWED = [
  "salmonella_enterica"
]

process seqsero {
    container "${LOCAL_REGISTRY}/bioinfo/seqsero2:1.1.1--e1654ae3eb"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 5.GB, task.attempt ) }
    cpus { taskCpus(32, task.attempt) }   
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'seqSero2/SeqSero_result.tsv', saveAs: { "${base}_result.tsv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: 'seqSero2/Extracted_antigen_alleles.fasta', saveAs: { "${base}_extracted_antigen_alleles.fasta" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: 'seqSero2/SeqSero_log.txt', saveAs: { "${base}_log.txt" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    when:
      isSpeciesSupported(genus_species, GENUS_SPECIES_ALLOWED, reads, task.process)    
    input:
      tuple val(riscd_input), path(reads)
      val(genus_species)
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_seqsero"
      is_fasta = r1.getName() ==~ /.+\.fa(sta)?$/
      if (is_fasta) {
        """
          SeqSero2_package.py -d seqSero2 -p ${task.cpus} -t 4 -i ${r1} -n ${md.cmp} -m k
        """     
      } else if (isIlluminaPaired(reads)) { 
        """    
          SeqSero2_package.py -d seqSero2 -p ${task.cpus} -t 2 -i ${r1} ${r2} -n ${md.cmp}
        """
      } else if (isIonTorrent(reads)) {
        """
          SeqSero2_package.py -d seqSero2 -p ${task.cpus} -t 3 -i ${r1} -n ${md.cmp}
        """                      
      }
}

workflow step_4TY_serotype__seqsero {
    take: 
      reads
      genus_species
    main:
      seqsero(reads, genus_species)
}

workflow {
    step_4TY_serotype__seqsero(getSingleInput(), param('genus_species'))
}