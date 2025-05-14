nextflow.enable.dsl=2

include { extractDsRef;getEmpty;flattenPath; parseMetadataFromFileName; executionMetadata; extractKey;taskMemory } from '../functions/common.nf'
include { getSingleInput;getReferences;getReferenceCodes;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent;isSegmentedMapping } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '2AS_mapping'
def METHOD = 'ivar' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def getExDt(reference, ex) {
    try {        
      reflist = getReferenceCodes()
      if (reflist.size() < 2) {
        return ex.dt
      } 
      if (!reference) {
        log.warn "${reference} should not be empty"
        return ex.dt
      }      
      if (!reflist.contains(reference)) {
        log.warn "${reference} not found in ${reflist}"
        return ex.dt
      }
      return ex.dt + reflist.indexOf(reference)
    } catch(Throwable t) {
        log.warn "getExDt: unexpected exception: ${t.asString()}"
        return ex.dt
    }
}

def canBeAggregated(actual) {
    try {        
      expected = getReferenceCodes()
      if (expected.size() < 2) {
          return false
      }
      if (actual.size() != expected.size()) {
        log.warn "expected ${expected.size()} references, got: ${actual.size()}"
        return false
      } 
      if (actual.any {!expected.contains(it)}) {
        log.warn "there is some unexpected value in ${actual} (values should be: ${expected})"
        return false
      }
      return true
    } catch(Throwable t) {
        log.warn "canBeAggregated: unexpected exception: ${t.asString()}"
        return false
    }
}

process snippy {
    container "ghcr.io/genpat-it/snippy:4.5.1--7be4a1c45a"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 8.GB, task.attempt ) }
    maxForks 4
    when:
      reference_path && reference_path.exists() && !reference_path.empty() && isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(reference_path)
    output:
      path "**/${base_ref}*"
      tuple path("snippy/${base_ref}.bam"), val(reference), emit: bam  
      path '*.sh', hidden: true
      path '*_input.json'
    afterScript "echo '${stepInputs([riscd_input,riscd_ref], md, [dt: ex_dt], STEP, METHOD, [reference:reference])}' > ${base_ref}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/result", pattern: '{**/*.tab,**/*.fa,**/*.bam*,**/*.bed,**/*.csv,**/*.vcf*,**/*.gff,**/*.html,**/*.txt}', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '**/*.log', saveAs: { filename -> flattenPath(filename) }     
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_snippy.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: "*.json"
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_vdsnippy_${reference}"
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

process samtools_pileup {
    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      tuple path(bam), val(reference)
    output:
      tuple path("${base_ref}.pileup"), val(reference), emit: pileup    
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools_pileup.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${METHOD}_${reference}"
      """
        #input should be DS10561267-DT200702_2020.TE.89540.1.2_vdsnippy_NC045512.bam
        samtools mpileup -d 1000 -A -Q 0 ${bam} > ${base_ref}.pileup
      """
}

process ivar {
    container "quay.io/biocontainers/ivar:1.3--h089eab3_1"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      tuple path(pileup), val(reference)
    output:
      path '*'
      path '*.sh', hidden: true
      tuple val(riscd), val(reference), path("${base_ref}.fasta"), emit: consensus    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/result", pattern: '*.fasta'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_${METHOD}.cfg" }
    script:
      md = parseMetadataFromFileName(pileup.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${METHOD}_${reference}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
        # Note : samtools mpileup output must be piped into ivar consensus
        cat ${pileup} | ivar consensus -p ${base_ref} -q ${params.step_2AS_mapping__ivar___q} -m 1 > ${base_ref}.log
        mv ${base_ref}.fa ${base_ref}.fasta
      """
}

process samtools_depth {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
        input:
      tuple path(bam), val(reference)
      val(method)
    output:
      tuple path("${base_ref}.coverage"), val(reference), emit: coverage
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools_depth.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${method}_${reference}"
      """
      samtools depth -a ${bam} | awk '{ if (\$3!=0) c++;s+=\$3}{h++} END { if (c!=0) print s/c; else print 0;if (h!=0) print c/h; else print 0 }' > ${base_ref}.coverage
	    """
}

process coverage_minmax {
    container "ghcr.io/genpat-it/samtools:0.1.19--f3869562fe"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 4.GB, task.attempt ) }
        input:
      tuple path(bam), val(reference)
      val(method)
    output:
      path '*.csv'
      path '*.sh', hidden: true
      tuple path("${base_ref}_samtools_depth.txt"), val(reference), emit: coverage_depth  
      tuple val(md.ds), val(reference), path("${base_ref}_import_coverage_minmax.csv"), emit: coverage_extra
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '{*.csv,*.txt}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_minmax.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${method}_${reference}"
      """
      samtools view -F 4 -c ${bam} > samtools_view.txt
      samtools depth ${bam} > ${base_ref}_samtools_depth.txt
      /scripts/coverage_minmax.py ${md.cmp} ${md.ds} samtools_view.txt ${base_ref}_samtools_depth.txt ${base_ref}_import_coverage_minmax.csv
	    """
}

process coverage_check {
    container "ghcr.io/genpat-it/python3:3.10.1--29cf21c1f1"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      tuple path(coverage), path(consensus), val(reference)
      val(context)
    output:
      tuple val(md.ds), val(reference), path("${base_ref}_import_coverage.csv"), emit: coverage_basic
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '{*.csv,*.check}'      
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_check.cfg" }
    script:
      //TODO fix output folder
      md = parseMetadataFromFileName(consensus.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_${context}_${reference}"
      coverages= coverage.toRealPath().toFile().readLines()
      """
      /scripts/coverage.py ${coverage} ${consensus} ${reference} ${md.cmp} ${md.ds.substring(2)} ${base_ref}.check ${base_ref}_import_coverage.csv 
	    """
}

process coverage_check_group {
  container "ubuntu:20.04"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory { taskMemory( 200.MB, task.attempt ) }
  when:
    (files instanceof java.util.Collection) && files.size() > 1  
  input:
    tuple val(key), path(files)
    val(method)
  output:
    path '*'
    path '*.sh', hidden: true
  publishDir mode: 'rellink', "${params.outdir}/${res_folder}/meta", pattern: '.command.sh', saveAs: { "${base}_${method}_coverage_check.cfg" }
  publishDir mode: 'rellink', "${params.outdir}/${res_folder}/result", pattern: '*.csv'
  script:
    def coverage_file = files.flatten()[0]
    md = parseMetadataFromFileName(coverage_file.getName())
    base = "${md.ds}-${ex.dt}_${md.cmp}"
    res_folder = isSegmentedMapping() ? "${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}" : "aggregate"
    referenceOrdered = getReferenceCodes()
    files.sort{ c1, c2 ->
      def ref1 = c1.getName().replaceAll(".+_${method}_", "") -~ /_import.+/ 
      def ref2 = c2.getName().replaceAll(".+_${method}_", "") -~ /_import.+/ 
      referenceOrdered.indexOf(ref1) <=> referenceOrdered.indexOf(ref2) 
    }
    """
    cat ${files} > ${base}_${method}_coverage_full.csv.tmp
    head -n 1 ${base}_${method}_coverage_full.csv.tmp > ${base}_${method}_coverage_full.csv
    grep -E "^[A-Z]" -v ${base}_${method}_coverage_full.csv.tmp >> ${base}_${method}_coverage_full.csv
    """  
}

process coverage_check_merge {
  container "ubuntu:20.04"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory { taskMemory( 200.MB, task.attempt ) }
  input:
    tuple val(key), val(reference), path(covMinMax), path(covBasic)
    val(method)
  output:
    tuple val(md.ds), path("${base_ref}_import_coverage_merged.csv"), emit: coverage_merged
    path '*'
    path '*.sh', hidden: true
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '*.csv'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_check_merge.cfg" }
  script:
    md = parseMetadataFromFileName(covMinMax.getName())
    ex_dt = getExDt(reference, ex)
    base = "${md.ds}-${ex_dt}_${md.cmp}"
    base_ref = "${base}_${method}_${reference}"
    """
    paste -d, ${covMinMax} ${covBasic} | cut -d, -f1,2,3,4,5,8,9,10,11,12,13 > ${base_ref}_import_coverage_merged.csv
    """  
}

process coverage_plot {
    container "quay.io/biocontainers/matplotlib:3.1.2"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    input:
      tuple path(coverage_depth), val(reference)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/result", pattern: '*.png'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex_dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_plot.cfg" }
    script:
      md = parseMetadataFromFileName(coverage_depth.getName())
      ex_dt = getExDt(reference, ex)
      base = "${md.ds}-${ex_dt}_${md.cmp}"
      base_ref = "${base}_ivar_${reference}"
      """
      /scripts/coverage_plot.py ${coverage_depth} ${base_ref}_coverage_plot.png
	    """
}

process aggregate {
  container "ubuntu:20.04"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory { taskMemory( 200.MB, task.attempt ) }
  when:
    isSegmentedMapping() && canBeAggregated(references)
  input:
    tuple val(riscd_input), val(references), path(consensus)
  output:
    path '*'
    path '*.sh', hidden: true
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_${METHOD}_aggregate.cfg" }
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: "*.json"
  afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [references:referenceOrdered])}' > ${base}_input.json"
  script:
    md = parseMetadataFromFileName(consensus[0].getName())
    base = "${md.ds}-${ex.dt}_${md.cmp}"
    sorted = []
    referenceOrdered = getReferenceCodes()
    references.eachWithIndex{ ref, index ->
        sorted[referenceOrdered.indexOf(ref)] = consensus[index]
    }
    """
      for f in ${sorted.join(" ")} ; do 
        if [ "0" -ne "`grep -cve '^\\s*\$' -e ">" \${f}`" ] ;
            then
                cat "\${f}" | awk 1 >> ${base}_${METHOD}_aggregate.fasta
        fi
      done
    """  
}

workflow step_2AS_mapping__ivar {
  take: 
    reads
    reference 
  main:
    bam = snippy(reads, reference).bam
    pileup = samtools_pileup(bam).pileup
    consensus = ivar(pileup).consensus

    coverage_minmax(bam, 'vdsnippy')
    coverage_minmax.out.coverage_depth | coverage_plot

    coverage = samtools_depth(bam, 'vdsnippy').coverage

    coverage.cross(consensus) { extractDsRef(it) }.map { 
        return [ it[0][0], it[1][2], it[0][1] ]
    }.set { coverageRefAndConsensus }
    coverageBasic = coverage_check(coverageRefAndConsensus, 'ivar').coverage_basic

    crossedChecks = coverage_minmax.out.coverage_extra.cross(coverageBasic) { it[0] + "-" + it[1] }
    .map { [ it[0][0], it[0][1], it[0][2], it[1][2] ] }

    coverage_check_group(coverage_check_merge(crossedChecks, 'vdsnippy').coverage_merged | groupTuple, 'vdsnippy')
    
    reads.cross(consensus | groupTuple) {{ extractKey(it) }}
      .map { [
         it[0][0], // riscd
         it[1][1], // reference codes
         it[1][2]  // consensus file
      ]}.set { consensus_to_aggregate}

    aggregate(consensus_to_aggregate)
  
  emit:
    consensus = consensus.map { it[0,2] }
    coverage_depth = coverage_minmax.out.coverage_depth
}

workflow {
    getSingleInput().cross(getReferences('any')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }

    step_2AS_mapping__ivar(input.reads, input.refs)
}