nextflow.enable.dsl=2

include { extractDsRef;parseMetadataFromFileName;executionMetadata;extractKey;taskMemory;getEmpty } from '../functions/common.nf'
include { isRunningFromSampleSheet } from '../functions/samplesheet.nf'
include { getInput;getReference;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'


def ex = executionMetadata()

def STEP = '2AS_indexing'
def METHOD = 'bowtie' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"


// Indexing with Bowtie2. TAKES: 1 path to fasta file; OUTPUTS: variable number of index files from Bowtie2 (*.bt2).
process bowtie_build {
	container "${LOCAL_REGISTRY}/bioinfo/bowtie2:2.1.0--37ad014737"
	input:
		path(fasta)
	output:
		path '**'
		path '{*.sh,*.log}', hidden: true
		afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
		publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta'
		publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.bt2'
		publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
		publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
		publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
	script:
		md = parseMetadataFromFileName(fasta.getName())
		base = "${md.ds}-${ex.dt}_${md.cmp}_bowtie_build"
		//reference command: bowtie2-build -f ${fasta} ${accession}
		"""
		bowtie2-build -f ${fasta} ${md.cmp}
		"""
}


// Workflow: TAKES from params "cmp" and "riscd" and RUNS process for `bowtie2-build`.
workflow step_2AS_indexing__bowtie {
	take:
		input_data
	main:
		fasta = input_data.map { it[1] }
		bowtie_build(fasta)
}


// Main Workflow: runs named workflow with `getInput()`'s nested tuple output as parameter from params file (necessary: NCBI accession number or genpat ID, riscd).
workflow {
	step_2AS_indexing__bowtie(getInput())
}