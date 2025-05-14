nextflow.enable.dsl=2

include { parseMetadataFromFileName;getRisCd;executionMetadata;taskMemory;extractKey;stepInputs;parseRISCD } from '../functions/common.nf'
include { getInputOf;getReferences;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent;param;optionalOrDefault } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '1PP_filtering'
def METHOD = 'samtools' 

process extractHeader {
    container 'quay.io/biocontainers/samtools:0.1.19--hf89b575_7'
    memory { taskMemory( 18.GB, task.attempt ) }
    input:
      tuple val(bam), path(bampath)
    output:
      stdout
    script:
      """
      samtools view -H ${bampath} | grep '^@SQ' | head -n 1 | cut -f2 | sed 's/SN://' | tr -d '\\n'
      """
}


process bedfile {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    memory { taskMemory( 250.MB, task.attempt ) }
    input:
      val(header)
    output:
      path("regions.bed"), emit: bed    
    script:
      regions = param('step_1PP_filtering__samtools__regions')     
      """
#!/usr/bin/env python3

try:
  n = open("regions.bed","w")
  regions = "${regions}"
  header = "${header}"
  for pos in regions.split(","):
    print(pos)
    start = pos.split("-")[0].strip().replace(".","")
    end = pos.split("-")[1].strip().replace(".","")
    n.write(f"{header}\t{start}\t{end}\\n")
  n.close()
except Exception as e:
  print(f"Error: {e}")
	    """  
}

process samtools {
    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 18.GB, task.attempt ) }
    input:
      path(bed)
      tuple val(riscd_input), path(bampath)
    output:
      tuple val(riscd), path('*.fastq.gz'), emit: filtered
      path '*'
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [reference:reference,filtered_region: regions])}' > ${base}_samtools_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "${base}_samtools_R*.fastq.gz"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}_samtools.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_samtools.cfg" }
    script:
      regions = param('step_1PP_filtering__samtools__regions')
      md = parseMetadataFromFileName(bampath.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      filteredR1 = "${base}_samtools_R1.fastq.gz"
      filteredR2 = "${base}_samtools_R2.fastq.gz"
      riscd = getRisCd(md, ex, STEP, METHOD)
      reference = (bampath.getName() =~ /([^_]+)\.bam$/)[0][1]
      type = optionalOrDefault('step_1PP_filtering__samtools__type', 'paired')
      if (type == 'paired') {
        """
        samtools view -b -L regions.bed ${bampath} | samtools fastq -1 ${filteredR1} -2 ${filteredR2} -0 /dev/null -s /dev/null -n
        """
}     else if (type == 'single') {
        """
        samtools view -b -L regions.bed ${bampath} | samtools fastq - > ${filteredR1}
        """
      }      
}

workflow step_1PP_filtering__samtools {
    take: 
      bam
    main:
      bed = bedfile(extractHeader(bam))
      samtools(bed, bam)
    emit:
      samtools.out.filtered
}

workflow {
    bam = getInputOf("*.bam")
    step_1PP_filtering__samtools(bam)
}

