nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;getGenusSpeciesOptionalUnkeyed } from '../functions/parameters.nf'
include { stepInputs;flattenPath;csv2map } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '4TY_lineage'
def METHOD = 'snapperdb' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def SNAPPER_CONFIGS = [
  brucella_melitensis : "brucella_melitensis.txt",
  brucella_abortus : "brucella_abortus.txt"
]

def GASTROSNAPPER_CONFPATH='/bioinfonas/databases/PROGRAMS/snapperdb/user_configs'
def GASTROSNAPPER_REFPATH='/bioinfonas/databases/PROGRAMS/snapperdb/reference_genomes'

def LOCKFILE = '/biowork/tmp/.SNAPPERDB_LOCKFILE'

process fastq_to_db {
    container "${LOCAL_REGISTRY}/bioinfo/snapper-db:1.0.6--bddce8087e"
    containerOptions = " -v ${GASTROSNAPPER_CONFPATH}:/user_configs:ro -v ${GASTROSNAPPER_REFPATH}:/reference_genomes:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    when:
      SNAPPER_CONFIGS.containsKey(genus_species?.toLowerCase())    
    input:
      tuple val(riscd_input), path(reads)
      val(genus_species)
      path(lockfile)
    output:
      tuple val(config), path("${base}.log"), emit: log
      val(genus_species), emit: genus_species
      path '**'
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [config: config])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*/*.vcf*', saveAs: { filename -> flattenPath(filename) }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (t1,t2) = reads
      config = SNAPPER_CONFIGS.get(genus_species?.toLowerCase())
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      sample_id = md.cmp.replaceAll(/\./, '-')
      """
        #!/bin/bash -euo pipefail -l
        run_snapperdb.py get_strains -c ${config} > strains.tsv
        CURRENT_ADDRESS=`grep "^${sample_id}\t" strains.tsv | wc -l || true`
        if [ "\${CURRENT_ADDRESS}" -eq "0" ]
        then
          mv ${t1} ${sample_id}.fastq.gz && mv ${t2} ${sample_id}_2.fastq.gz
          run_snapperdb.py fastq_to_vcf -g . -c ${config} ${sample_id}.fastq.gz ${sample_id}_2.fastq.gz &> ${base}.log
          ls snpdb/*.filtered.vcf
          flock ${lockfile} -c 'run_snapperdb.py vcf_to_db -c ${config} snpdb/*.filtered.vcf' &>> ${base}.log
        else
          echo "[${md.cmp}] Sample '${sample_id}' already added for config: ${config}" >> ${base}.log    
        fi
      """
}

process update {
  container "${LOCAL_REGISTRY}/bioinfo/snapper-db:1.0.6--bddce8087e"
  containerOptions = " -v ${GASTROSNAPPER_CONFPATH}:/user_configs:ro -v ${GASTROSNAPPER_REFPATH}:/reference_genomes:ro"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory params.max_memory
  when:
    !params.step_4TY_lineage__snapperdb___skip_update
  input:
    tuple val(config), path(logfile)
    path(lockfile)
  output:
    tuple val(config), path("${base}.tsv"), path("strains.tsv"), emit: snp_address, optional: true
    path '*'
    path '*.sh', hidden: true
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.tsv'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.log'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
  script:
    md = parseMetadataFromFileName(logfile.getName())
    base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}_update"    
    sample_id = md.cmp.replaceAll(/\./, '-')
    result = ''
    for (line in logfile.toRealPath().toFile().readLines()) {
      def matcher = (line =~ /.* depth is ([\d\.]+) - cut+off is ([\d\.]+)/)
      if (matcher.matches()) {
        def depth = matcher.group(1) as double
        def cutoff = matcher.group(2) as double
        if (depth < cutoff) {
          log.warn "Sample '${md.cmp}' not added to: ${config}. Depth: ${depth}, cutoff: ${cutoff}"
        } else {
          result = 'JUST_ADDED'
        }
        break
      } else if (line ==~ /.* Sample '[^']+' already added .*/) {
        log.warn "Sample '${md.cmp}' already added for config: ${config}"
        result = 'ALREADY_ADDED'
        break
      } else if (line ==~ /.* is already in SNPdb strains_snps .*/) {
        log.warn "Sample '${md.cmp}' is already in SNPdb strains_snps, but no SNP address was assigned; config: ${config}"
        break
      }
    } 
    if (result == 'JUST_ADDED') {
      """
        #!/bin/bash -euo pipefail -l
        flock ${lockfile} -c 'run_snapperdb.py update_distance_matrix -c ${config}' &> ${base}.log
        flock ${lockfile} -c 'run_snapperdb.py get_strains -c ${config}' > strains.tsv
        echo -e "genpat\tsnp_address" > ${base}.tsv
        OUTLIER=`grep "Outlier Z-Score" ${base}.log | wc -l || true`
        if [ "\${OUTLIER}" -eq "0" ]
        then
          grep -P '^${md.cmp}\t' strains.tsv | head -n 1 >> ${base}.tsv
        else
          echo -e "${sample_id}\tTBD" >> ${base}.tsv
        fi
      """ 
    } else if (result == 'ALREADY_ADDED') {
      """
        #!/bin/bash -euo pipefail -l
        flock ${lockfile} -c 'run_snapperdb.py get_strains -c ${config}' > strains.tsv
        echo -e "genpat\tsnp_address" > ${base}.tsv
        grep -P '^${md.cmp}\t' strains.tsv | head -n 1 >> ${base}.tsv
      """ 
    } else {
      def message = "snp_address not calculated"
      log.warn "${message}"
      """
        #!/bin/bash -euo pipefail -l
        echo -e "genpat\tsnp_address" > ${base}.tsv
        echo -e "${sample_id}\tNA" >> ${base}.tsv
        echo "${message}" > ${base}.log
      """ 
    }
}

process outbreak_check {
  container "ubuntu:20.04"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory { taskMemory( 500.MB, task.attempt ) }
  input:
    tuple val(config), path(snp_address), path(snp_addresses)
  output:
    path '*.sh', hidden: true
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.log'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
  script:
    md = parseMetadataFromFileName(snp_address.getName())
    base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}_check"    

    fullAddress = csv2map(snp_address, "\\t").snp_address
    sixthLevelAddress = fullAddress -~ /[^\.]+$/

    sample_id = md.cmp.replaceAll(/\./, '-')

    level7group = snp_addresses.toRealPath().toFile().findAll { it.contains('\t') && it.split('\t')[1] == fullAddress}.collect { it.split('\t')[0] } - sample_id
    level6group = snp_addresses.toRealPath().toFile().findAll { it.contains('\t') && it.split('\t')[1].startsWith(sixthLevelAddress)}.collect { it.split('\t')[0] } - sample_id

    (level6message, level7message) = [ '', '' ]
    
    if (level7group.size()) {
      level7message =  "[${md.cmp}] has the same full SNP address [${fullAddress}] of: ${level7group}"
      log.warn level7message
    }
    if (level6group.size()) {
      level6message = "[${md.cmp}] has the same SNP address prefix [${sixthLevelAddress}*] of: ${level6group}"
      log.warn level6message
    }
    """
      [ -z "${level7message}" ] || echo "${level7message}" >> ${base}.log 
      [ -z "${level6message}" ] || echo "${level6message}" >> ${base}.log 
    """
}

workflow step_4TY_lineage__snapperdb {
  take: 
    reads
    genus_species      
  main:
    result = fastq_to_db(reads.collect(), genus_species, LOCKFILE).log
    update(result, LOCKFILE).snp_address | outbreak_check
}

workflow {
  step_4TY_lineage__snapperdb(getSingleInput(), getGenusSpeciesOptionalUnkeyed())
}

workflow update_distance_matrix {
  update(getGenusSpeciesOptionalUnkeyed(), LOCKFILE)
}
