nextflow.enable.dsl=2

//1 import functions from functions/common.nf and functions/parameters.nf
include { stepInputs;taskMemory;parseMetadataFromFileName;executionMetadata } from '../functions/common.nf'
include { getInputOf;parseRISCD;param } from '../functions/parameters.nf'

//2 define variables for construction of output directory structure and output file names:
//2.1 define 'ex' - execution date metadata: a dictionary of key value pairs (example [anno:2024, dt:DT240906])
def ex = executionMetadata()

//2.2 define parts of the analysis name used to name output files and subdirectories
def STEP = '4TY_distance'
def METHOD = 'mash_dist' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

//3 define variables for name and location of reference files (optional, only for fixed reference).
//  Reference needs to be located in ~/ngsmanager/assets/${ENTRYPOINT}
def referenceCode = 'dataset_sketch'
def referencePath = "${params.assets_dir}/${ENTRYPOINT}/dataset_sketch.msh"
def referenceRiscd = '240903-020240903125047730-4TY_distance-mashdist'

//4 process definition
process mash_dist {
    container "${LOCAL_REGISTRY}/bioinfo/mash:2.3--debdd7eb23"
    memory { taskMemory(8.GB, task.attempt) }
    cpus 8
    input:
        //4.1 takes a tuple from input channel and gets paths for the 2 elements of the tuple to use as input
        tuple path(sample), path(reference)
    output:
        //4.2 rules for ouptut file locations
        //4.3.1 sets location for logs and sh files at root and sets them as hidden (adds '.' to file name)
        path '*'
        path '{*.sh,*.log}', hidden: true
        //4.3.2 runs bash code after analysis to output an input report. echo redirects to a json the output of
        //      the stepInputs() function (a dictionary with process and input info)
        afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
        //4.3.3 creates a directory structure for the output inside the 'results' directory (${params.outdir}).
        //      Files produced by the workflow will be symlinked (mode: 'rellink') as specified. 'pattern'
        //      defines which files will be matched.
        publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.tsv'
        publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.txt'
        publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*_input.json'
        publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
        publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
        //4.4 runs the script for each sample/iteration (element passed from channel)
        //4.4.1 defines the 'md' variable ('metadata'): a dictionary of substrings from the input sample's name.
        //      It is used to build the base output file suffix ('base' variable) and to build the directory
        //      structure in the output block.
        md = parseMetadataFromFileName(sample.getName())
        base = "${md.ds}-${ex.dt}_${md.cmp}_mash_dist"
        //4.4.2 Actual script to run. If there is the need to use bash variables, the '$' must be escaped ('\${var}').
        """
        mash dist ${reference} ${sample} > ${base}.tsv
        """
//        filename="${sample}"
//        IFS='_' read -r part1 part2 part3 <<< "\$filename"
//        echo "mash dist ${reference} ${sample} > \${part2}_mash_dist.tsv" >> testecho${base}.txt
}

//5 named workflow
workflow step_4TY_distance__mashdist {
    //5.1 definition of input: the item passed as input to the process is named and massed to the 'main' block
    take:
        sketch
    //5.2 processes the input item before passing it to the channel(s) to feed the process call
    main:
        //5.2.1 declaration of new variable to store the reference's path (working directory)
        reference = referencePath
        //5.2.2 for each item in input (samples), it extracts the path to the file to to use as input ('sample'),
        //      then for each of them it creates a new variable (a tuple), pairing each sample with the reference.
        //      The tuples are set to allow usage outside the indented block.
        sketch.map { entry -> 
            def sample = entry[1][0]
            tuple(sample,reference)
        }.set { sampleWithReference }
        //sampleWithReference.view()
        //5.2.3 feeds the channel with sample-reference pairs to the defined process
        mash_dist(sampleWithReference)
}

//6 default unnamed workflow. It's always executed first if it's present.
//  We only use it to call the named workflow, passing as input the output of the getInput() function
//  (or its variants). The output of getInput() is a nested list, containing the data string of the analysis
//  that produced the files used as input and the list of files produced by it (inner list).
workflow {
    step_4TY_distance__mashdist(getInputOf('*.msh'))
}
