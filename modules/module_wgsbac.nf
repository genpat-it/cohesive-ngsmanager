nextflow.enable.dsl=2

include { executionMetadata; taskMemory; extractKey; isSpeciesSupported; parseMetadataFromFileName; taskCpus } from '../functions/common'
include { stepInputs; getRisCd } from '../functions/common.nf'
include { param; getInput; getReference; isCompatibleWithSeqType } from '../functions/parameters.nf'

def db_kraken = "/bioinfonas/databases/PROGRAMS/kraken/minikraken_8GB_20200312"
def db_confindr = "/bioinfonas/databases/PROGRAMS/confindr/confindr_full"

def GENUS_SPECIES_ALLOWED = ["Salmonella", "Brucella"]


// named process 1 - creates input tsv for process 2:
//                 * expects 2 IN channels, and 1 element from each channel
//                 * emits 1 OUT channel, returning 1 file
//                 * performs tuple unpacking of 2 elements of list received from reads channel
process table_builder {
	container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
	memory { taskMemory( 250.MB, task.attempt ) }
	input:
		path(reads)
		path(ref_path)
	output:
		path("metadata.tsv"), emit: tsv
	script:
		(r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
		"""
		#!/usr/bin/env python

		# Input data
		## Take input reads from channel
		sample = "${r1}".split('_R1')[0]
		reads = [[sample, "${r1}", "${r2}"]]
		
		## Take reference data from channel
		ref_path = "${ref_path}"
		ref = ref_path.split('/')[-1]

		# Define output file
		output_file = 'metadata.tsv'
		summary_file = 'summary.tsv'

		# Create the table structure
		table = [["ID", "Common.name", "Type", "File_left_reads", "File_right_reads", "File_assembly", "Lat.Long.of.Isolation"]]

		# Add samples
		for e in reads:
			sample_id = e[0].split('-')[0]
			left = sample_id + '_R1' + '.fastq.gz'
			right = sample_id + '_R2' + '.fastq.gz'
			table.append([sample_id, sample_id, "new", left, right, "", ""])

		# Add reference genome
		if ref:
			## Extract ID from filename
			ref_id = ref.split('.')[0].split('-')[0]
			ref_fasta = ref_id + '.fasta'
			table.append([ref_id, ref_id, "reference", "", "", ref_fasta, ""])
		else:
			raise ValueError("reference file is required!")

		# Write to TSV file
		with open(output_file, "w") as out:
			for i in table:
				print(*i, sep="\t", file=out)
		
		## Write summary
		with open(summary_file, 'w') as smry:
			print("${r1}", left, sep="\t", file=smry)
			print("${r2}", right, sep="\t", file=smry)
			print(ref, ref_fasta, sep="\t", file=smry)
		"""
}


// named process 2 - wgsbac docker container
//      * options: mounts database directories
//      * only executed for Illumina paired-end && for supported species
//      * expects 4 IN channels, and 1 element from each channel
//         * ref_path + genus_species needed to link input files in process dir
process iaea_wgsbac {
	container "${LOCAL_REGISTRY}/bioinfo/wgsbac:1.0--007111424d"
	containerOptions "-u 0:0 -v ${db_confindr}:/mnt/confindr -v ${db_kraken}:/mnt/minikraken_8GB_20200312 -v ${task.workDir}:/wd"
	memory { taskMemory( 24.GB, task.attempt ) }
	cpus params.max_cpus
	stageInMode 'copy'
//	when:
//		isCompatibleWithSeqType(rawreads, ['illumina_paired'], task.process) && isSpeciesSupported(genus_species, GENUS_SPECIES_ALLOWED, rawreads, task.process)
	input:
		path(tsv)
		path(reads)
		path(ref_path)
		val genus
	output:
		path '*'
		path '*.sh', hidden: true
		//publishDir mode: 'rellink', "wgsbac_results/", pattern: "wgsbac_results/*"
		publishDir mode: 'rellink', "result/wgsbac_results/", pattern: 'wgsbac_results/*'
		publishDir mode: 'rellink', "meta", pattern: '*_input.json'

	script:
		//md = parseMetadataFromFileName(fq_r1)
		//base = "${md.cmp}_module_IAEA-WGSBAC"

		if (genus == 'Salmonella') {
			"""
			#!/usr/bin/bash

			set -euo pipefail

			mkdir -p wgsbac_results/input
			mkdir -p wgsbac_results/finalAssembly
			for i in DS* ; do
				if [[ \$i =~ DS.+_R1.+ ]]; then
					mv \$i  wgsbac_results/input/\${i/-*/}_R1.fastq.gz
				elif [[ \$i =~ DS.+_R2.+ ]]; then
					mv  \$i wgsbac_results/input/\${i/-*/}_R2.fastq.gz
				elif [[ \$i =~ DS.+.fasta ]]; then
					mv  \$i wgsbac_results/finalAssembly/\${i/-*/}.fasta
				fi
			done

			wgsbac.pl --table ${tsv} --results wgsbac_results --kraken /mnt/minikraken_8GB_20200312 --mlst senterica --seqsero --sistr --amr --amrorganism Salmonella --virulence --spi --conf /mnt/confindr --snippy --plasmid --cpus ${task.cpus} --run
			"""
		} else if (genus == 'Brucella') {
			"""
			#!/usr/bin/bash

			set -euo pipefail

			mkdir -p wgsbac_results/input
			mkdir -p wgsbac_results/finalAssembly
			for i in DS* ; do
				if [[ \$i =~ DS.+_R1.+ ]]; then
					mv \$i  wgsbac_results/input/\${i/-*/}_R1.fastq.gz
				elif [[ \$i =~ DS.+_R2.+ ]]; then
					mv  \$i wgsbac_results/input/\${i/-*/}_R2.fastq.gz
				elif [[ \$i =~ DS.+.fasta ]]; then
					mv  \$i wgsbac_results/finalAssembly/\${i/-*/}.fasta
				fi
			done

			wgsbac.pl --table ${tsv} --results wgsbac_results --kraken /mnt/minikraken_8GB_20200312 --mlst brucella --virulence --conf /mnt/confindr --snippy --plasmid --cpus ${task.cpus} --run
			"""
		} else {
			error "Invalid genus provided: only Salmonella and Brucella are allowed!" 
		}
}


// named workflow - takes 3 channels and:
//                1) extracts list of paths to R1 and R2 + path to reference
//                2) runs first and second processes passing 2 channels each
workflow pipeline_wgsbac {
	take:
		rawreads
		genus_species
		reference
	main:
		reads = rawreads.map { it[1] }
		//genus_species_ch = Channel.value(genus_species)
		ref_path = reference.map { it[3] }

		def (genus, species) = genus_species.contains("_") ? genus_species.split('_') : [ genus_species, null ]

		tables = table_builder(reads, ref_path)
		iaea_wgsbac(tables[0], reads, ref_path, genus)

		//summary = vcf_input.collectFile { [ "summary.tsv", it[0] + '\t' + filename(it[1]) + '\n' ] }
}


// main workflow - runs pipeline_wgsbac workflow by passing as channels:
//               * list of input data ([cmp, [/path/to/R1, /path/to/R2]])
//               * genus of MO ('genus_species' key from params.json)
//               * list of reference data ([DS, cmp, AN, /path/to/fasta])
workflow {
	pipeline_wgsbac(getInput(), param('genus_species'), getReference('fa'))
}
