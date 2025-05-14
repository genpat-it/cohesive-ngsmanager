nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;param } from '../functions/parameters.nf'

def ex = executionMetadata()

def STEP = '4AN_AMR'
def METHOD = 'filtering'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process abricate_filtering {
    container "${LOCAL_REGISTRY}/bioinfo/python3:3.10.1--29cf21c1f1"
    memory { taskMemory( 250.MB, task.attempt ) }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(calls)
      val(coverage)
      val(identity)
    output:
      path '*'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [identity: identity, coverage: coverage])}' > ${base}_abricate_filtering_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'    
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_abricate_filtering.cfg" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*_abricate_calls.txt'    
    script:
      md = parseMetadataFromFileName(calls.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """    
        #!/usr/bin/env python3

        new = open("${base}_COV${coverage}_ID${identity}_abricate_calls.txt", "w")
        res = open("${calls}", 'r').readlines()
        new.write(res[0])
        for l in res:
            try:
                if float(l.split("\t")[9]) >= ${coverage} and float(l.split("\t")[10]) >= ${identity}:
                    new.write(l.replace(l.split("\t")[0],"${calls}"))
            except:
                pass
        new.close()
      """
}

workflow step_4AN_AMR__filtering {
    take: 
      data
      coverage
      identity
    main:
      abricate_filtering(data, coverage, identity);
}

workflow {
    step_4AN_AMR__filtering(getSingleInput(), param('coverage'), param('identity'))
}
