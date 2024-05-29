nextflow.enable.dsl=2

include { logHeader;taskMemory;taskTime } from '../functions/common.nf'
include { getNCBICodes } from '../functions/parameters.nf'

log.info logHeader('NGSMANAGER')

def REF_FORMATS = ['fasta', 'gb']

def ASSEMBLY_FASTA_SUFFIX='_genomic.fna.gz'
def ASSEMBLY_GB_SUFFIX='_genomic.gbff.gz'
def ASSEMBLY_MAX_DOWNLOAD_SIZE= params.ncbi_max_download_assembly
def SRA_MAX_DOWNLOAD_SIZE= params.ncbi_max_download_sra

process ncbi_assembly_refseq {
    container "ncbi/edirect:12.5"
    tag "${ref}"
    memory { taskMemory( 250.MB, task.attempt ) }
    time { taskTime( 10.m, task.attempt ) }        
    maxForks 3
    when:
      ref ==~ /^GCF.+$/ 
    input:
      val(ref)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${ref}/fasta/result", pattern: '{*.fasta}'    
    publishDir mode: 'rellink', "${params.outdir}/${ref}/fasta/meta", pattern: '{*.json,*.xml,*.txt,*.log}'
    publishDir mode: 'rellink', "${params.outdir}/${ref}/fasta/meta", pattern: '.command.sh', saveAs: { "${ref}.cfg" }
    script:
      assert ref ==~ /^GCF_[\w\.]+$/
      """
        esearch -db assembly -query "${ref} AND latest_refseq [PROP]" -sort Accession | esummary > search_refseq_summary       
        DOCID=`xtract -input search_refseq_summary -pattern DocumentSummary -position last -first Id | awk 'NF'`
        VERSION=`xtract -input search_refseq_summary -pattern DocumentSummary -position last -first Synonym/RefSeq | awk 'NF'`
        if  [[ ! -z "\${DOCID}" ]] && [[ ! -z "\${VERSION}" ]] ;
        then
          ANAME=`xtract -input search_refseq_summary -pattern DocumentSummary -position last -first AssemblyName | awk 'NF'`
          echo -e "\${ANAME}\n\${VERSION}" > \${VERSION}.txt          
          efetch -db assembly -id \${DOCID} -format docsum > \${VERSION}.xml
          efetch -db assembly -id \${DOCID} -format docsum -mode json > \${VERSION}.json
          FTP_PATH=`xtract -input \${VERSION}.xml -pattern DocumentSummary -position last -first FtpPath_RefSeq | awk 'NF'`
          FTP_BASENAME=`echo \${FTP_PATH} | sed 's/.*\\///g'`
          echo "\${FTP_PATH}/\${FTP_BASENAME}${ASSEMBLY_FASTA_SUFFIX}" > \${VERSION}.log
          curl --max-filesize ${ASSEMBLY_MAX_DOWNLOAD_SIZE} \${FTP_PATH}/\${FTP_BASENAME}${ASSEMBLY_FASTA_SUFFIX} | gunzip -c > \${VERSION}.fasta
        else 
          (>&2 echo "${ref}: not found or not the latest version"; exit 2) 
        fi      
      """
}

process ncbi_assembly_genbank {
    container "ncbi/edirect:12.5"
    tag "${ref}"
    memory { taskMemory( 250.MB, task.attempt ) }
    time { taskTime( 10.m, task.attempt ) }         
    maxForks 3
    when:
      ref ==~ /^GCA.+$/     
    input:
      val(ref)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${ref}/gb/result", pattern: '{*.gb}'      
    publishDir mode: 'rellink', "${params.outdir}/${ref}/gb/meta", pattern: '{*.json,*.xml,*.txt,*.log}'
    publishDir mode: 'rellink', "${params.outdir}/${ref}/gb/meta", pattern: '.command.sh', saveAs: { "${ref}.cfg" }
    script:
      assert ref ==~ /^GCA_[\w\.]+$/
      """
        esearch -db assembly -query "${ref} AND latest_genbank [PROP]" -sort Accession | esummary > search_genbank_summary       
        DOCID=`xtract -input search_genbank_summary -pattern DocumentSummary -position last -first Id | awk 'NF'`
        VERSION=`xtract -input search_genbank_summary -pattern DocumentSummary -position last -first Synonym/Genbank | awk 'NF'`
        if  [[ ! -z "\${DOCID}" ]] && [[ ! -z "\${VERSION}" ]] ;
        then
          ANAME=`xtract -input search_genbank_summary -pattern DocumentSummary -position last -first AssemblyName | awk 'NF'`
          echo -e "\${ANAME}\n\${VERSION}" > \${VERSION}.txt            
          efetch -db assembly -id \${DOCID} -format docsum > \${VERSION}.xml
          efetch -db assembly -id \${DOCID} -format docsum -mode json > \${VERSION}.json
          FTP_PATH=`xtract -input \${VERSION}.xml -pattern DocumentSummary -position last -first FtpPath_GenBank | awk 'NF'`
          FTP_BASENAME=`echo \${FTP_PATH} | sed 's/.*\\///g'`
          echo "\${FTP_PATH}/\${FTP_BASENAME}${ASSEMBLY_GB_SUFFIX}" > \${VERSION}.log
          curl --max-filesize ${ASSEMBLY_MAX_DOWNLOAD_SIZE} \${FTP_PATH}/\${FTP_BASENAME}${ASSEMBLY_GB_SUFFIX} | gunzip -c > \${VERSION}.gb
        else 
          (>&2 echo "${ref}: not found or not the latest version"; exit 2) 
        fi    
      """
}

process ncbi_sra {
    container "quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    tag "${code}"
    memory { taskMemory( 2.GB, task.attempt ) }
    time { taskTime( 20.m, task.attempt ) }         
    maxForks 3
    input:
      val(code)
    output:
      path '*'
      path '{*.sh,*.log}', hidden: true
    publishDir  mode: 'rellink', "${params.outdir}/${code}/fastq/result", pattern: '*.gz'
    publishDir  mode: 'rellink', "${params.outdir}/${code}/fastq/meta", pattern: '{*.json,*.txt}'
    publishDir  mode: 'rellink', "${params.outdir}/${code}/fastq/meta", pattern: '.command.log', saveAs: { "${code}.log" }
    publishDir  mode: 'rellink', "${params.outdir}/${code}/fastq/meta", pattern: '.command.sh', saveAs: { "${code}.cfg" }
    script:
      """
        prefetch "${code}" -H 1 -X ${SRA_MAX_DOWNLOAD_SIZE} |& tee .size_check.tmp
        [ "`grep -c 'is larger than maximum allowed' .size_check.tmp`" == "0" ]
        vdb-dump "${code}" --info -f json >  "${code}.json"
        fasterq-dump -L info --skip-technical --split-files "${code}"      
        for f in *.fastq; do gzip \$f; done 
        echo "${code}\n${code}" >  "${code}.txt"
      """
}


process ncbi_nuccore {
    container "ncbi/edirect:12.5"
    tag "${ref}-${format}"
    memory { taskMemory( 250.MB, task.attempt ) }
    time { taskTime( 5.m, task.attempt ) }         
    maxForks 3
    input:
      val(ref)
      each format 
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir  mode: 'rellink', "${params.outdir}/${ref}/${format}/result", pattern: '{*.fa*,*.gb}'
    publishDir  mode: 'rellink', "${params.outdir}/${ref}/${format}/meta", pattern: '{*.json,*.txt}'
    publishDir mode: 'rellink', "${params.outdir}/${ref}/${format}/meta", pattern: '.command.sh', saveAs: { "${ref}.cfg" }
    script:
      efetch_format = format == 'gb' ? 'gbwithparts' : format
      """
        efetch -db nuccore -id ${ref} -format docsum -mode json > ${ref}.json
        efetch -db nuccore -id ${ref} -format url > ${ref}.txt
        efetch -db nuccore -id ${ref} -format acc >> ${ref}.txt
        [[ `cat ${ref}.txt | wc -l` -eq '2' ]] || (>&2 echo "Accession Number or URL not retrieved!"; exit 1)    
        efetch -db nuccore -id ${ref} -format ${efetch_format} > ${ref}.${format} && grep '\\S' ${ref}.${format}
      """
}

workflow {
  getNCBICodes().branch {
      gcf: it ==~ /(?i)^GCF.+/
      gca: it ==~ /(?i)^GCA.+/
      sra: it ==~ /(?i)^SRR.+|(?i)^ERR.+|(?i)^DRR.+/
      nuccore: true
  }.set { branched }
  ncbi_assembly_refseq(branched.gcf)
  ncbi_assembly_genbank(branched.gca)
  ncbi_nuccore(branched.nuccore, REF_FORMATS)
  ncbi_sra(branched.sra)
}