nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata; extractKey;taskMemory;stepInputs } from '../functions/common.nf'
include { optionalBoolean;getSingleInput;getReference;getReferenceCodes;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '2AS_mapping'
def METHOD = 'snippy' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def PUBLISH_SAMTOOLS_DEPTH = optionalBoolean('step_2AS_mapping__snippy__publish_samtools_depth')

process snippy {
    container "ghcr.io/genpat-it/snippy:4.5.1--7be4a1c45a"
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
      tuple path("snippy/${base_ref}.bam"), val(reference), emit: bam
      tuple path("snippy/${base_ref}.vcf"), val(reference), emit: vcf
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
      base_ref = "${base}_${METHOD}_${reference}"
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

process samtools {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      tuple path(bam), val(reference)
    output:
      tuple path("${base_ref}_view.txt"), path("${base_ref}_depth.txt"), val(reference), emit: mapping_data
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "${base_ref}_depth.txt", enabled: PUBLISH_SAMTOOLS_DEPTH
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_samtools_${reference}"
      """
      samtools view -F 4 -c ${bam} > ${base_ref}_view.txt
      samtools depth -a ${bam} > ${base_ref}_depth.txt
	    """
}

process coverage {
    container "ghcr.io/genpat-it/samtools:0.1.19--f3869562fe"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 4.GB, task.attempt ) }
    input:
      tuple path(view), path(depth), val(reference)
    output:
      path '*.csv'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      md = parseMetadataFromFileName(view.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_coverage_${reference}"
      """
      awk '{ if (\$3!=0) c++;s+=\$3}{h++} END { if (c!=0) print s/c; else print 0;if (h!=0) print c/h; else print 0 }' ${depth} > coverage.txt
      /scripts/coverage_minmax.py ${md.cmp} ${md.ds} ${view} ${depth} ${base_ref}_import_coverage_minmax.csv
      /scripts/coverage_light.py coverage.txt ${reference} ${md.cmp} ${md.ds.substring(2)} ${base_ref}_import_coverage.csv 
      paste -d, ${base_ref}_import_coverage_minmax.csv ${base_ref}_import_coverage.csv  | cut -d, -f1,2,3,4,5,8,9,10,11,12,13 > ${base_ref}_import_coverage_merged.csv      
      """
}

process coverage_plot {
    container "quay.io/biocontainers/matplotlib:3.1.2"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    input:
      tuple path(view), path(depth), val(reference)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.png'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base_ref}.log" }    
    script:
      md = parseMetadataFromFileName(depth.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_pyplot_${reference}"
      """
      /scripts/coverage_plot.py ${depth} ${base_ref}_coverage_plot.png
	    """
}

process bcftools {
    container "quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      tuple path(vcf), val(reference)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*.vcf,*.txt,*.tsv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base_ref}.log" }
    script:
      md = parseMetadataFromFileName(vcf.getName())
            base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_bcftools_${reference}"
      """
      bcftools norm -a ${vcf} -o decomposed_complex.vcf.tmp
      bcftools view -v snps decomposed_complex.vcf.tmp -o snps_view.vcf.tmp
      bcftools norm -m +any snps_view.vcf.tmp -o ${base_ref}_snps_only.vcf
      bcftools view -v indels decomposed_complex.vcf.tmp -o ${base_ref}_indels_only.vcf
      bcftools stats ${vcf} | grep '^SN' | sed 's/number of \\(.*\\):/\\1/g' | awk -v FS='\t' -v OFS='\t' '{print \$3"\t"\$4}' > ${base_ref}_stats_summary.tsv
      bcftools stats ${base_ref}_snps_only.vcf | grep '^SN' | sed 's/number of \\(.*\\):/\\1/g' | awk -v FS='\t' -v OFS='\t' '{print \$3"\t"\$4}' > ${base_ref}_snps_only_stats_summary.tsv
      grep -v '^#' ${base_ref}_snps_only.vcf | wc -l > ${base_ref}_snps_count.txt
	    """
}

workflow step_2AS_mapping__snippy {
  take: 
    reads
    reference 
  main:
    bam = snippy(reads, reference).bam
    mapping_data = samtools(bam).mapping_data 
    coverage_plot(mapping_data)
    coverage(mapping_data)
    bcftools(snippy.out.vcf)
}

workflow {
    getSingleInput().cross(getReference('any')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }
    step_2AS_mapping__snippy(input.reads, input.refs)
}