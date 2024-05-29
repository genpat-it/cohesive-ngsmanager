nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata; extractKey;taskMemory;stepInputs } from '../functions/common.nf'
include { getSingleInput;getReference;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '2AS_mapping'
def METHOD = 'snippy' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process snippy {
    container "${LOCAL_REGISTRY}/bioinfo/snippy:4.5.1--7be4a1c45a"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 8.GB, task.attempt ) }
    maxForks 4
    when:
      reference_path && reference_path.exists() && !reference_path.empty() && (isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process))
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(reference_path)
    output:
      path "**/${base_ref}*"
      path '*_input.json'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_ref], md, ex, STEP, METHOD, [reference:reference, ref_format:ref_format])}' > ${base_ref}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{**/*.tab,**/*.fa,**/*.bam*,**/*.bed,**/*.csv,**/*.vcf*,**/*.gff,**/*.html,**/*.txt}', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: "*.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base_ref}.log" }    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_snippy_${reference}"
      ref_format = (reference_path.getName() ==~ /.+\.f(n)?(a)?(sta)?$/) ? 'fasta' : 'gb'
      is_fasta = r1.getName() ==~ /.+\.fa(sta)?$/
      if (is_fasta) {
        """
        trap "find -name "*.sam" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --ctgs ${r1} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        NUMVAR=\$((`cat snippy/*.tab | wc -l`-1))
        echo -e "Sample\tRef\tNumVar" >> ${base_ref}.check
        echo -e "${md.cmp}\t${reference}\t\$NUMVAR" >> ${base_ref}.check
        """     
      } else if (isIlluminaPaired(reads)) {
        """
        trap "find -name "*.sam" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --R1 ${r1} --R2 ${r2} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        NUMVAR=\$((`cat snippy/*.tab | wc -l`-1))
        echo -e "Sample\tRef\tNumVar" >> ${base_ref}.check
        echo -e "${md.cmp}\t${reference}\t\$NUMVAR" >> ${base_ref}.check
        """
      } else if (isIonTorrent(reads)) {
        """
        trap "find -name "*.sam" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --se ${r1} --outdir snippy --prefix ${base_ref} --quiet  --tmpdir tmp  &> ${base_ref}-cli.log
        NUMVAR=\$((`cat snippy/*.tab | wc -l`-1))
        echo -e "Sample\tRef\tNumVar" >> ${base_ref}.check
        echo -e "${md.cmp}\t${reference}\t\$NUMVAR" >> ${base_ref}.check
        """      
      }
}

workflow step_2AS_mapping__snippy {
  take: 
    reads
    reference 
  main:
    snippy(reads, reference)
}

workflow {
    getSingleInput().cross(getReference('any')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }
    step_2AS_mapping__snippy(input.reads, input.refs)
}