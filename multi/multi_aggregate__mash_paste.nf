nextflow.enable.dsl=2

include { taskMemory } from '../functions/common.nf'
include { getInputOf;param } from '../functions/parameters.nf'


process mash_paste {
	container "${LOCAL_REGISTRY}/bioinfo/mash:2.3--debdd7eb23"
	memory { taskMemory( 8.GB, task.attempt ) }
	cpus 8
	input:
		path(sketches)
	output:
		path '*'
		path 'dataset_sketch.msh', emit: refsketch
		path '{*.sh,*.log}', hidden: true
		publishDir mode: 'rellink', "${params.outdir}/result", pattern: '{dataset_sketch.msh}'
		publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '.command.sh', saveAs: { "mash_paste.cfg" }
		publishDir mode: 'rellink', "${params.outdir}/meta", pattern: '{dataset_sketch_info.txt}'
	script:
		"""
			mash paste dataset_sketch ${sketches}
			mash info dataset_sketch.msh > dataset_sketch_info.txt
		"""
}

workflow multi_aggregate__mash {
	take:
		input
	main:
		sketches = input.collect { it[1] }
		mash_paste(sketches)
		//refsketch = mash_paste(sketches).refsketch
		//pairs = refsketch.combine(input)
		//test(pairs)
		//mash_dist(pairs)
		//pairs.view()
}

workflow {
	multi_aggregate__mash(getInputOf('*.msh'))
}