nextflow.enable.dsl=2

include { step_4TY_plasmid__mobsuite } from '../steps/step_4TY_plasmid__mobsuite'
include { isSpeciesSupported;taskMemory;flattenPath;parseMetadataFromFileName;executionMetadata } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent;param } from '../functions/parameters.nf'

POINFINDER_GENUS_ALLOWED = [
  'campylobacter',
  'enterococcus_faecalis',
  'enterococcus_faecium',
  'escherichia_coli',
  'helicobacter_pylori',
  'klebsiella',
  'mycobacterium_tuberculosis',
  'neisseria_gonorrhoeae',
  'plasmodium_falciparum',
  'salmonella',
  'staphylococcus_aureus'
]

def ex = executionMetadata()

def getPointfinderParam(gsp) {
 try {  
    def genus_species = gsp ? gsp.toLowerCase() : ''
    def (genus, species) = genus_species.contains("_") ? genus_species.split('_') : [ genus_species, null ]
    if (POINFINDER_GENUS_ALLOWED.contains(genus_species)) {
        // genus_species allowed
        return genus_species;
    }       
    if (POINFINDER_GENUS_ALLOWED.contains(genus)) {
        // ALL genus allowed
        return genus
    }
    return ''
  } catch(Throwable t) {
      exit 1, "unexpected exception: ${t.asString()}"
  } 
}

process staramr {
    container "ghcr.io/genpat-it/staramr:0.9.1--8fe6b5a239"
    containerOptions = "--user root"
    tag "${riscd_input}/${plasmid_name}"
    memory { taskMemory( 2.GB, task.attempt ) }
    input:
      tuple val(riscd_input), path(plasmid_file)
      val(genus_species)
    output:
      path '**'
      path '{*.sh,*.log}', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/plasmids/${plasmid_name}/result", pattern: 'result/{*.tsv,*.xlsx}', saveAs: { filename -> "${base}_${flattenPath(filename)}" }  
    publishDir mode: 'rellink', "${params.outdir}/plasmids/${plasmid_name}/result", pattern: 'result/hits/*', saveAs: { filename -> "hits/${base}_hits_${flattenPath(filename) -~ /_DS.+/}.fasta" }  
    publishDir mode: 'rellink', "${params.outdir}/plasmids/${plasmid_name}/meta", pattern: 'result/*.txt', saveAs: { filename -> "${base}_${flattenPath(filename)}" }  
    publishDir mode: 'rellink', "${params.outdir}/plasmids/${plasmid_name}/meta", pattern: '*.json'
    publishDir mode: 'rellink', "${params.outdir}/plasmids/${plasmid_name}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
    publishDir mode: 'rellink', "${params.outdir}/plasmids/${plasmid_name}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      plasmid_name = plasmid_file.getName() -~ /\.fa.*$/
      base = "${plasmid_name}_staramr"
      pointfinder_db = getPointfinderParam(genus_species)
      extra = pointfinder_db ? " --pointfinder-organism ${pointfinder_db}" : ''
      """
        sed 's/^>.*/>/g' ${plasmid_file} | awk '{for(i=1;i<=NF;i++){if(\$i~/^>/){\$i=">contig"++count}}} 1' > ${base}.fasta
        staramr search ${extra} -o result ${base}.fasta
      """
}

workflow module_plasmids {
    take: 
        assembly
        genus_species
    main:
        plasmids = step_4TY_plasmid__mobsuite(assembly).plasmids

        plasmids.multiMap { it ->
            riscd: it[0]
            plasmids: it[1]
        }.set { branched }

        input = branched.riscd.combine(branched.plasmids.flatten())
        staramr(input, genus_species)
}

workflow {
    module_plasmids(getSingleInput(), param('genus_species'))
}