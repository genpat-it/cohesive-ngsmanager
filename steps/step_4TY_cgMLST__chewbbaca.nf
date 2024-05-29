nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;taskTime } from '../functions/common.nf'
include { getInput;param;optionalOrDefault;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

SPECIES_SCHEMA = [
  listeria_monocytogenes : ['l_mono_chewie_1748_220623'],
  escherichia_coli : ['e_coli_chewie_2360_210531'],
  salmonella_enterica : ['s_enterica_chewie_3255_210531']
]

SCHEMAS = [
  l_mono_chewie_1748_220623 : "/schemas/Listeria_monocytogenes_Pasteur_cgMLST_2022-06-23T18_03_54.613576.zip",
  e_coli_chewie_2360_210531 : "/schemas/Escherichia_coli_INNUENDO_wgMLST_2021-05-31T14_24_05.304225.zip",
  s_enterica_chewie_3255_210531 : "/schemas/Salmonella_enterica_INNUENDO_cgMLST_2021-05-31T20_28_21.350919.zip"
]

CHEWBBACA_SINGLE_END_PARAMS = [
  'l_mono_chewie_1748_220623': ' --minimum-length 144 --st 0.1 --bsr 0.6 ',
  'e_coli_chewie_2360_210531': ' --minimum-length 0 --st 0.01 --bsr 0.6 --genes-list /schemas/Escherichia_coli_INNUENDO_cgMLST_EFSA_filterlist.txt ',
  's_enterica_chewie_3255_210531': ' --minimum-length 0 --st 0.01 --bsr 0.6 --genes-list /schemas/Salmonella_enterica_INNUENDO_cgMLST_EFSA_filterlist.txt '
]

CHEWBBACA_PAIRED_END_PARAMS = [
  'l_mono_chewie_1748_220623': ' --minimum-length 144 ',
  'e_coli_chewie_2360_210531': ' --minimum-length 0 --genes-list /schemas/Escherichia_coli_INNUENDO_cgMLST_EFSA_filterlist.txt ',
  's_enterica_chewie_3255_210531': ' --minimum-length 0 --genes-list /schemas/Salmonella_enterica_INNUENDO_cgMLST_EFSA_filterlist.txt '
]

def getExtraParams(schema) {
  try {   
    if (isIonTorrent(null)) {
      if (CHEWBBACA_SINGLE_END_PARAMS.containsKey(schema)) {
        return CHEWBBACA_SINGLE_END_PARAMS[schema]
      }
    } else {
      if (CHEWBBACA_PAIRED_END_PARAMS.containsKey(schema)) {
        return CHEWBBACA_PAIRED_END_PARAMS[schema]
      }
    }        
    return ''
  } catch(Throwable t) {
      exit 1, "getExtraParams, unexpected exception: ${t.asString()}"
  } 
}

def ex = executionMetadata()

def STEP = '4TY_cgMLST'
def METHOD = 'chewbbaca' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def getSchema(gsp, schema) {
  try {  
    def genus_species = gsp ? gsp.toLowerCase() : ''
    def (genus, species) = genus_species.contains("_") ? genus_species.split('_') : [ genus_species, null ]

    def allowedSchemas = []
    if (SPECIES_SCHEMA.containsKey(genus_species)) {
      allowedSchemas =  SPECIES_SCHEMA.get(genus_species) 
    } else if (SPECIES_SCHEMA.containsKey(genus)) {
      allowedSchemas = SPECIES_SCHEMA.get(genus)
    }
    if (!allowedSchemas) {
      log.debug "no compatible schemas found for genus_species: ${genus_species}"
      return null
    }
    if (!schema) {
      return allowedSchemas[0] //at this stage, use the first allowed schema
    }
    if (allowedSchemas.contains(schema)) {
      return schema
    }
    log.warn "schema ${schema} not compatible with genus_species: ${genus_species}"
    return null;      
  } catch(Throwable t) {
      exit 1, "unexpected exception: ${t.asString()}"
  } 
}

process chewbbaca {
    container "${LOCAL_REGISTRY}/bioinfo/chewbbaca-w-chewie-schemas:2.8.5--16b816c96d"
    memory { taskMemory( 8.GB, task.attempt ) }
    time { taskTime( 30.m, task.attempt ) }    
    cpus { [32, params.max_cpus as int].min() }
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      getSchema(genus_species, schema)
    input:
      tuple val(riscd_input), path(assembly)
      val genus_species
      val schema
    output:
      path '**'
      path("${base}_results_statistics.tsv"), emit: stats
      tuple path("${base}_results_alleles.tsv"), path('schema/'), emit: alleles
      tuple path("${base}_results_alleles.tsv"), path("${base}_new_alleles.txt"), val(schemaName), emit: alleles_with_new
      path '{*.sh,*.log}', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [schema:schemaName, chewbbaca: '2.8.5'])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.tsv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "${base}_new_alleles.txt"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: 'listGenes2Call.txt', saveAs: { "${base}_listGenes2Call.txt" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      md = parseMetadataFromFileName(assembly.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      schemaName = getSchema(genus_species, schema)
      schemaPath = SCHEMAS.get(schemaName) 
      newAlleleKey = assembly.getName().replaceAll('_', '-')
      speciesSpecificParams = getExtraParams(schemaName)
      """
        #!/bin/bash -euo pipefail
        unzip ${schemaPath} -d schema > /dev/null
        chmod -R 777 schema
        mkdir input && cp ${assembly} input/
        chewBBACA.py AlleleCall -i input -g schema -o results --cpu ${task.cpus} --force-continue --verbose ${speciesSpecificParams}
        grep "${newAlleleKey}" schema/*.fasta -A1 -h | grep -v "\\-\\-" > ${base}_new_alleles.txt || echo "no INF alleles found"
        mv results/*/results_alleles.tsv ${base}_results_alleles.tsv
        mv results/*/results_contigsInfo.tsv ${base}_results_contigsInfo.tsv
        mv results/*/results_statistics.tsv ${base}_results_statistics.tsv        
      """
}

process chewbbaca_check {
    container "quay.io/biocontainers/python:3.9"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"    
    memory { taskMemory( 1.GB, task.attempt ) }
    time { taskTime( 5.m, task.attempt ) }    
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      path(chewbbacaStats)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
      path("${base}_import_chewbbaca_check.csv"), emit: check
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.log', saveAs: { "${base}_chewbbaca_check.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_chewbbaca_check.cfg" }

    script:
      md = parseMetadataFromFileName(chewbbacaStats.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      /scripts/chewieCheck.py --stat ${chewbbacaStats} > ${base}_import_chewbbaca_check.csv      
      """ 
}

process hashing {
    container "${LOCAL_REGISTRY}/bioinfo/hashing:1.0--29180a232f"    
    memory { taskMemory( 2.GB, task.attempt ) }
    time { taskTime( 10.m, task.attempt ) }    
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple path(chewbbaca_result), path(schema_path)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '_hashed_results.tsv', saveAs: { "${base}_chewbbaca_results_crc32.tsv" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}_hashing.log" }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_hashing.cfg" }
    script:
      md = parseMetadataFromFileName(chewbbaca_result.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
        mask_matrix.py -i ${chewbbaca_result} -o masked_results.tsv
        alleleprofile_hasher.py -p masked_results.tsv -d ${schema_path} -o ./_hashed_results.tsv   
        rm -Rf ${schema_path}/* 
      """
}

workflow step_4TY_cgMLST__chewbbaca {
    take: 
      assembly
      genus_species
      schema
    main:
      //assume channels are already crossed
      chewbbaca_result = chewbbaca(assembly, genus_species, schema)      
      hashing(chewbbaca_result.alleles)
      chewbbaca_check(chewbbaca_result.stats).check
}

workflow {
    step_4TY_cgMLST__chewbbaca(getInput(), param('genus_species'), optionalOrDefault('schema', ''))
}
