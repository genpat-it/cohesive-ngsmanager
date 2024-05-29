nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;_getSingleReference } from '../functions/parameters.nf'
include { stepInputs;csv2map;getRisCd;parseRISCD } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_lineage'
def METHOD = 'westnile' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

STEP_ASSETS_DIR = "${params.assets_dir}/${ENTRYPOINT}"

def HCOV_THRESHOLD = params.step_4TY_lineage__westnile___threshold

def WESTNILE_LINEAGE_REFERENCES_PATH = "${STEP_ASSETS_DIR}/westnile_lineage_references.json"

WESTNILE_LINEAGE_REFERENCES = new groovy.json.JsonSlurper().parseText(file(WESTNILE_LINEAGE_REFERENCES_PATH).text)


def getReferenceForLineage(lineageFile) {
   try {        
      def lineage = csv2map(lineageFile, ",").lineage
      assert lineage
      def refData = WESTNILE_LINEAGE_REFERENCES.find { it.lineage == lineage }
      assert refData
      return [ refData.ref_riscd, refData.ref_code, "${STEP_ASSETS_DIR}/${refData.ref_path}" ]
   } catch(Throwable t) {
       exit 1, "unexpected exception: ${t.asString()}"
   }   
}

def isValidLineage(lineageResult) {
   try {        
      def lineage = csv2map(lineageResult[1], ",").lineage   
      def md = parseMetadataFromFileName(lineageResult[1].getName())
      def notAssigned =  lineage == 'NA'
      def notDetermined = lineage == 'ND'
      if (notAssigned) {
        log.warn "[${md.cmp}] Lineage not assigned"
      } else if (notDetermined) {
        log.warn "[${md.cmp}] Lineage not determined"
      }
      return !notAssigned && !notDetermined
   } catch(Throwable t) {
       exit 1, "unexpected exception: ${t.asString()}"
   }   
}
 
def getWNVReferences() {
   try {        
      references = Channel.fromList(WESTNILE_LINEAGE_REFERENCES.collect { [ it.ref_riscd, it.ref_code, "${STEP_ASSETS_DIR}/${it.ref_path}" ] } )
   } catch(Throwable t) {
      exit 1, "unexpected exception: ${t.asString()}"
   }   
}

process snippy {
    container "${LOCAL_REGISTRY}/bioinfo/snippy:4.5.1--7be4a1c45a"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 8.GB, task.attempt ) }
    maxForks 4
    when:
      reference_path && reference_path.exists() && !reference_path.empty()
    input:
      tuple val(riscd_input), path(reads), val(riscd_ref), val(reference), path(reference_path)
    output:
      path "**/${base_ref}*"
      path '*'
      tuple val(riscd), val(reference), path("snippy/${base_ref}.bam"),  emit: bam  
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${METHOD}_${reference}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      if (r2) {
        """
        trap "find -name "*.sam" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --R1 ${r1} --R2 ${r2} --outdir snippy --prefix ${base_ref} --quiet --tmpdir tmp  &> ${base_ref}.log
        
        """
      } else {
        """
        trap "find -name "*.sam" -delete ; rm -Rf tmp" EXIT
        mkdir tmp
        snippy --reference ${reference_path} --ctgs ${r1} --outdir snippy --prefix ${base_ref} --quiet --tmpdir tmp  &> ${base_ref}.log
        """
      }
}

process samtools_depth {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 500.MB, task.attempt ) }
    input:
      tuple val(riscd_input), val(reference), path(bam)
    output:
      path '*'
      tuple val(riscd), path("${base_ref}.coverage"), emit: coverage
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      lineage = WESTNILE_LINEAGE_REFERENCES.find { it.ref_code == reference }.lineage
      assert lineage
      base_ref = "${base}_${METHOD}_${reference}"
      riscd = getRisCd(md, ex, STEP, METHOD)
      """
      samtools depth -a ${bam} | awk '{ if (\$3!=0) c++;s+=\$3}{h++} END { if (c!=0) print s/c; else print 0;if (h!=0) print c/h; else print 0 }' > ${base_ref}.coverage
	    echo "${reference}" >> ${base_ref}.coverage
	    echo "${lineage}" >> ${base_ref}.coverage
      """
}

process lineage_selection {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 250.MB, task.attempt ) }
    input:
      tuple val(riscd_input), path(covfile)
    output:
      path '*'
      tuple val(riscd), path("${base}.csv"), emit: lineage, optional: true
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: 'errors.log', saveAs: { "${base}.err" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.csv'
    script:
      coverage_file = (covfile instanceof java.util.Collection) ? covfile.flatten()[0] : covfile
      md = parseMetadataFromFileName(coverage_file.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}_lineage"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
        for f in  ${covfile}; do cat \${f}  | tr '\n' ',' | cut -d ',' -f 2,3,4 ; done | sort -nr > ${base}_summary.csv
        /scripts/select_lineage.py --file_list ${covfile} --threshold ${HCOV_THRESHOLD} --output ${base}.csv &> ${base}.log
	    """
}

workflow step_4TY_lineage__westnile {
    take: 
      reads
    main:
      references = getWNVReferences()
      comb = reads.combine(references)
      bam = snippy(comb).bam
      coverages = samtools_depth(bam).coverage.groupTuple()
      lineage = lineage_selection(coverages).lineage 
    emit:
      lineage = lineage.filter { isValidLineage(it) }
}

workflow {  
    step_4TY_lineage__westnile(getSingleInput())
}