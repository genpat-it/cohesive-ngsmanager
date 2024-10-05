include { getEmpty;parseRISCD;flattenPath;executionMetadata;getGB;getRisCd } from './common.nf'

def getTrimmedReads(optional){
    def allowedAcc = ['1PP_trimming', '1PP_hostdepl', '1PP_generated']
    if (!params.containsKey('cmp') || !params.containsKey('trimming-riscd') || !(params['trimming-riscd'] instanceof CharSequence)) {
        if (optional) {
            log.debug "no trimmed reads available"
            return Channel.empty()
        }
        exit 2, "missing required params (cmp,trimming-riscd)";
    }
    def riscd =  params['trimming-riscd']
    def cmp = params.cmp.replaceAll("-", ".")
    def anno = cmp.substring(0,4)
    def md = parseRISCD(riscd)            
    if (!(md.acc in allowedAcc)) {
        exit 2, "unexpected acc value: ${md.acc}, expected: one of ${allowedAcc}"
    }      
    def path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/*.fastq*"
    def resChannel = Channel.fromPath(
            path, checkIfExists: params.check_file_existence && !optional
    )        
    .toSortedList( { a, b -> flattenPath(a.getName()) <=> flattenPath(b.getName()) } )             
    .collect()
    .map { [ riscd , it ] }
    resChannel.ifEmpty { 
        log.warn("file not found: '${path}'")
    }
    return resChannel
}

def getAssembly() {
    def allowedAcc = ['2AS_denovo', '2AS_import', '2AS_mapping']
    if (!params.containsKey('cmp') || !params.containsKey('riscd')) {
        exit 2, "missing required params (cmp,riscd)";
    }
    def cmp = params.cmp.replaceAll("-", ".")
    def anno = cmp.substring(0,4)
    def md = parseRISCD(params.riscd)       
    def path, resChannel     
    if (!(md.acc in allowedAcc)) {
        exit 2, "unexpected acc value: ${md.acc}, expected: one of ${allowedAcc}"
    }         
    if (md.acc == '2AS_import' || md.acc == '2AS_mapping') {
        path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/*.fasta"
        resChannel = Channel.fromPath(
                path, checkIfExists: params.check_file_existence
        )        
        .first()
        .collect()
        .map { [ params.riscd, it ] }

    } else {
        def pattern = (md.met == 'shovill' ? "*.fasta" : "*L?00.fasta")
        path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${pattern}"
        resChannel = Channel.fromPath(
                path, checkIfExists: params.check_file_existence
        )        
        .first()
        .collect()
        .map { [ params.riscd, it ] }
    }
    resChannel.ifEmpty { 
        log.warn("file not found: '${path}'")
    }
    return resChannel
}

def getDepletedReads(){
    def allowedAcc = ['1PP_hostdepl']
    if (!params.containsKey('cmp') || !params.containsKey('riscd')) {
        exit 2, "missing required params (cmp,riscd)";
    }
    def cmp = params.cmp.replaceAll("-", ".")
    def anno = cmp.substring(0,4)
    def md = parseRISCD(params.riscd)            
    if (!(md.acc in allowedAcc)) {
        log.error "unexpected acc value: ${md.acc}, expected: one of ${allowedAcc}"
        exit 2, "unexpected acc value: ${md.acc}, expected: one of ${allowedAcc}"
    }        
    def path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/*${md.ref ? '_' : ''}${md.ref}.fastq*"
    
    def resChannel = Channel.fromPath(
            path, checkIfExists: params.check_file_existence
    )        
    .toSortedList( { a, b -> flattenPath(a.getName()) <=> flattenPath(b.getName()) } )             
    .collect()
    .map { [ params.riscd, it ] }
    resChannel.ifEmpty { 
        log.warn("file not found: '${path}'")
    }
    return resChannel
}

def getReferenceCodes() {
    if (!params.containsKey('references') || !(params.references instanceof ArrayList)) {
        return []
    }
    return params.references.collect { it.code }
}

def getNCBICodes() {
    def val = param('reference')
    if (val instanceof ArrayList) {
        return Channel.fromList(val)
    } else {
        return Channel.fromList(val.tokenize(',\s\n\t\r'))
    }
}

def getReferences(type) {
    _getReferences(false, type, true, true)
}

def getReference(type) {
    _getReferences(false, type, true, false)
}

def getReferenceUnkeyed(type) {
    _getReferences(false, type, false, false).map { [ it[1], it[2], it[3] ]}
}

def getReferenceOptional(type) {
    _getReferences(true, type, true, false)
}

def _getReferences(optional, type, keyed, multi) {
    def crossValue = ""
    if (keyed) {
        // CROSSING PARAM
        if (!params.containsKey('riscd')) {
            exit 2, "one of: [ds, riscd] param should be provided";
        } 
        crossValue = parseRISCD(params.riscd).ds         
    }
    if (params.containsKey('references')) {
        assert params.references instanceof ArrayList
        if (!params.references && !optional) {
            exit 2, "No reference provided";
        }
        def refs = params.references.collect { refOriginal ->
            def ref = [:]
            refOriginal.keySet().each {
                ref["ref_${it}"] = refOriginal[it]
            }
            return ref
        }
        def result = refs.inject (Channel.empty()) {
            res, val -> res.concat(_getSingleReference(optional, type, val, crossValue))
        }              
        if (!multi && params.references.size() > 1) {
            log.warn "${params.references.size()} references provided; only the first one will be considered."
        }
        // never returns more than one reference if multi is false
        return multi ? result : result.first()
    } else {
        _getSingleReference(optional, type, params, crossValue)
    } 
}

def _getSingleReference(optional, type, refInput, crossValue) {
    def allowedAcc = ['2AS_mapping', '2AS_denovo', '2AS_import']
    if (refInput.containsKey('ref_cmp') && refInput.containsKey('ref_code') &&
          (
            (type == 'gb' && refInput.containsKey('ref_riscd_genbank') && refInput.ref_riscd_genbank)
              ||
            (type == 'fa' && refInput.containsKey('ref_riscd') &&  refInput.ref_riscd)
              ||
            (type == 'any' && (
                (refInput.containsKey('ref_riscd') &&  refInput.ref_riscd) || (refInput.containsKey('ref_riscd_genbank') && refInput.ref_riscd_genbank)
                              )
            )
          )
        ) {
        def useFastaRef = (type == 'fa') || (type == 'any' && (refInput.containsKey('ref_riscd') &&  refInput.ref_riscd))

        def cmp = refInput.ref_cmp.replaceAll("-", ".")
        def anno = cmp.substring(0,4)
        def filePattern, md, refRiscd 
        if (useFastaRef) {
            refRiscd = refInput.ref_riscd
            md = parseRISCD(refInput.ref_riscd)  
            filePattern = md.met == 'spades' ? '*L?00.fa*' : '*.fa*'
        } else {
            allowedAcc = ['4AN_genes']
            refRiscd = refInput.ref_riscd_genbank
            md = parseRISCD(refInput.ref_riscd_genbank)  
            filePattern = '*.gb*'
        }
        if (!(md.acc in allowedAcc)) {
            exit 2, "unexpected acc value: ${md.acc}, expected: one of ${allowedAcc}"
        }
        def path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${filePattern}"
        return Channel.fromPath(
             path, checkIfExists: params.check_file_existence && !optional
        )
        .first().map { [ crossValue, refRiscd, refInput.ref_code, it ] }
        .ifEmpty([ crossValue, '-', '-', getEmpty() ])       
    } else if (!optional) {
        exit 2, "could not find reference of type '${type}'";
    } else {
        log.debug "Not using (optional) reference, crossValue: ${crossValue}, type: ${type}"
        return Channel.of([ crossValue, '-', '-', getEmpty() ]) 
    }
}

def getHost() {
    _getParam('host', false, true)
}

def getHostOptional() {
    _getParam('host', true, true)
}

def getHostUnkeyed() {
    _getParam('host', false, false)
}

//https://www.gitmemory.com/issue/nextflow-io/nextflow/1388/564021238
def getGenusSpeciesOptionalUnkeyed(){
    _getParam('genus_species', true, false)
}

def getGenusSpeciesOptional(){
    _getParam('genus_species', true, true)
}

def getGenusSpecies(){
    _getParam('genus_species', false, true)
}

def getSpecies(){
    _getParam('species', false, true)
}

def getModule(){
    _getParam('module', false, true)
}

def getAbricateDatabase(){
    _getParam('abricate_database', false, true)
}

def getBlastDatabase(){
    _getParam('blast_database', false, true)
}

def getBlastDatabaseUnkeyed(){
    _getParam('blast_database', false, false)
}

def getKingdom(){
    _getParam('kingdom', false, true)
}

def getTaxIdsUnkeyed(){
    _getParam('taxids', true, false)
}

def getParamTaxaId() {
    _getParamAsValue('taxaid', false, null)
}

def getParamIncludeParents() {
    _getParamAsValue('include-parents', true, false)
}

def getParamIncludeChildren() {
    _getParamAsValue('include-children', true, false)
}

def _getParamAsValue(paramName, optional, defaultValue) {
    def res
    if (paramName && params.containsKey(paramName) && params[paramName] != '' && params[paramName] != null) {
        res = params[paramName]
    } else if (optional) {      
        res = defaultValue
    } else  {
        exit 2, "missing required param: $paramName";
    }  
    // no malicious chars inside
    if (res != null && res ==~ /(?s).*$[;|&><\(\)\n].*/) {
        exit 2, "possibile malicious content for param '${paramName}', value: '${value}'"
    }
    return res
}

def _getParam(paramName, optional, keyed){
     def crossValue
     if (keyed) {
        // CROSSING PARAM
        if (!params.containsKey('riscd')) {
            exit 2, "missing required param: riscd";
        }
        crossValue = parseRISCD(params.riscd).ds
        if (!crossValue) {
            exit 2, "could not get DS from riscd: ${params.riscd}"
        }
    }
    if (paramName && params.containsKey(paramName) && params[paramName]) {
        if (crossValue) {          
            return Channel.of( [ crossValue, params[paramName] ] )
        } else {
            return Channel.of( params[paramName] )
        }
    } else if (optional) {      
        if (crossValue) {       
            return Channel.of( [ crossValue, null ] )
        } else {
            return Channel.of( null )
        }
    } else  {
        exit 2, "missing required param: $paramName";
    }
}

def getDS() {
    if (!params.containsKey('riscd')) {
        exit 2, "missing required params (riscd)";
    }    
    return parseRISCD(params.riscd).ds
}

def isFullOutput() {
    return (params.containsKey('full_output') && params.full_output) || (params.containsKey('fullOutput') && params.fullOutput);
}

def getResult(cmp, riscd, filePattern) {
    def anno = cmp.substring(0,4)
    def md = parseRISCD(riscd)       
    def resChannel, path

    if (filePattern != null) {
        path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${filePattern}"
        resChannel = Channel.fromPath(
                path, checkIfExists: params.check_file_existence
        )        
        .collect()
        .map { [ riscd, it ] }        
    } else {
        if (md.acc in ['0SQ_rawreads', '1PP_hostdepl', '1PP_trimming', '1PP_filtering', '1PP_downsampling', '1PP_generated']) {
            def pattern = "*.fastq*"
            path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${pattern}"
            resChannel = Channel.fromPath(
                    path, checkIfExists: params.check_file_existence
            )        
            .toSortedList( { a, b -> flattenPath(a.getName()) <=> flattenPath(b.getName()) } )             
            .collect()
            .map { [ riscd, it ] }
        } else if (md.acc in ['2AS_import', '2AS_mapping']) {
            def pattern = "*.fasta"
            path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${pattern}"
            resChannel = Channel.fromPath(
                    path, checkIfExists: params.check_file_existence
            )        
            .first()
            .collect()
            .map { [ riscd, it ] }        
        } else if (md.acc in ['2AS_hybrid']) {
            def pattern = "*.fasta"
            path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${pattern}"
            resChannel = Channel.fromPath(
                    path, checkIfExists: params.check_file_existence
            )        
            .first()
            .collect()
            .map { [ riscd, it ] }        
        } else if (md.acc in ['2AS_denovo','2MG_denovo']) {
            def pattern = (md.met in ['spades','plasmidspades','unicycler','metaspades'] ? "*L?00.fasta" : "*.fasta")
            path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${pattern}"
            resChannel = Channel.fromPath(
                    path, checkIfExists: params.check_file_existence
            )        
            .first()
            .collect()
            .map { [ riscd, it ] }
        } else if (md.acc in ['4AN_genes']) {
            def pattern = "*.gff"
            path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${pattern}"
            resChannel = Channel.fromPath(
                    path, checkIfExists: params.check_file_existence
            )        
            .first()
            .collect()
            .map { [ riscd, it ] }        
        } else if (md.acc in ['4TY_cgMLST', '4TY_wgMLST'] && md.met in ['chewbbaca']) {
            def pattern = "*_results_${params.allelic_profile_encoding}.?sv"
            path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${pattern}"
            resChannel = Channel.fromPath(
                    path, checkIfExists: params.check_file_existence
            )        
            .first()
            .collect()
            .map { [ riscd, it ] }                
        } else if (md.acc in ['4AN_AMR']) {
            def pattern = "*abricate_calls.txt"
            path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${pattern}"
            resChannel = Channel.fromPath(
                    path, checkIfExists: params.check_file_existence
            )        
            .first()
            .collect()
            .map { [ riscd, it ] }        
        } else {
            exit 2, "unexpected acc value: ${md.acc}"
        }    
    }
    resChannel.ifEmpty { 
        log.warn("file not found: '${path}'")
    }
    return resChannel     
}

def _getAlleles(cmp, riscd, schema) {
    def md = parseRISCD(riscd)           
    def filePattern = schema ? "*_results_${schema}.?sv" : "*results_alleles.?sv"
    def anno = cmp.substring(0,4)
    def path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${filePattern}"
    def resChannel = Channel.fromPath(path, checkIfExists: params.check_file_existence)
    resChannel.ifEmpty { 
        log.warn("file not found: '${path}'")
    }
    return resChannel
}

def getAllelicProfile(cmp, riscd, schema) {
    def md = parseRISCD(riscd)           
    def filePattern = schema ? "*_results_${schema}.?sv" : "*results_alleles.?sv"
    def anno = cmp.substring(0,4)
    def path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/${filePattern}"
    def resChannel = Channel.fromPath(path, checkIfExists: params.check_file_existence)
    resChannel.ifEmpty { 
        log.warn("file not found: '${path}'")
    }
    return resChannel.map { [ riscd, it ] }
}

def getKrakenResults() {
    def allowedAcc = ['3TX_class']
    if (!params.containsKey('cmp') || !params.containsKey('riscd')) {
        exit 2, "missing required params (cmp,riscd)";
    }
    def md = parseRISCD(params.riscd)            
    if (!(md.acc in allowedAcc)) {
        exit 2, "unexpected acc value: ${md.acc}, expected: one of ${allowedAcc}"
    }     
    def cmp = params.cmp.replaceAll("-", ".")
    def anno = cmp.substring(0,4)    
    def path = "${params.inputdir}/${anno}/${cmp}/${md.acc}/${md.ds}-${md.dt}_${md.met}/result/*_kraken.{tsv,txt.gz}"
    def resChannel = Channel.fromPath(
            path, checkIfExists: params.check_file_existence
    )
    .collect()
    .map { [ params.riscd, it ] }
    resChannel.ifEmpty { 
        log.warn("file not found: '${path}'")
    }
    return resChannel
}

def getParam(paramName, optional, keyed) {
    _getParam(paramName, optional, keyed)
}

def getInputOf(pattern) {
    if (params.containsKey('input')) {
        assert params.input instanceof ArrayList
        params.input.inject (Channel.empty()) {
            res, val -> res.mix(getResult(val.cmp, val.riscd, pattern))
        }        
    } else if (params.containsKey('cmp') && params.containsKey('riscd')) {
        getResult(params.cmp, params.riscd, pattern)
    } else {
        exit 2, "missing required params (cmp,riscd) or (input)";
    } 
}

def getInput() {
    if (params.containsKey('input')) {
        assert params.input instanceof ArrayList
        params.input.inject (Channel.empty()) {
            res, val -> if (val.cmp) {
                res.mix(getResult(val.cmp, val.riscd, null))
            } else {
                res.mix(getInputFromPaths(val.key, val.paths))
            }
        }        
    } else if (params.containsKey('cmp') && params.containsKey('riscd')) {
        getResult(params.cmp, params.riscd, null)
    } else {
        exit 2, "missing required params (cmp,riscd) or (input)";
    } 
}

def getInputFolders() {
    getInputOf('')
}

def getVCFs() {
    getInputOf('{*[0-9].vcf,*.var.flt.vcf}')
}

def getSingleInput() {
    if (!params.containsKey('cmp') || !params.containsKey('riscd')) {
        exit 2, "missing required params (cmp,riscd)";
    }
    getResult(params.cmp, params.riscd, null)
}

def hasEnoughFastqData(rawreads) {
  try {
    def (r1,r2, r3, r4) = (rawreads instanceof java.util.Collection) ? rawreads : [rawreads,null]
    def sizeR1 = r1.toRealPath().countFastq()
    def sizeR2 = r2 ? r2.toRealPath().countFastq() : null
    log.debug "${r1}, size: ${sizeR1}"
    if (sizeR1 <= params.raw_reads_threshold && (sizeR2 == null || sizeR2 <= params.raw_reads_threshold)) {
      log.warn "Insufficient number of reads after trimming in: ${r1} (${sizeR1}) ${r2 ?: '-'} (${sizeR2 == null ? '-' : sizeR2}) (threshold: ${params.raw_reads_threshold})"
      return false;
    }
    return true;
  } catch(Throwable t) {
      exit 1, "unexpected exception: ${t.asString()}"
  } 
}

def hasFastqData(rawreads) {
  try {
    def (r1,r2) = (rawreads instanceof java.util.Collection) ? rawreads : [rawreads,null]
    def sizeR1 = r1.toRealPath().countFastq()
    def sizeR2 = r2 ? r2.toRealPath().countFastq() : null
    log.debug "${r1}, size: ${sizeR1}"
    if (sizeR1 == 0 && (sizeR2 == null || sizeR2 == 0)) {
      log.warn "No processable reads found in: ${r1} ${r2}"
      return false;
    }
    return true;
  } catch(Throwable t) {
      exit 1, "unexpected exception: ${t.asString()}"
  } 
}

def param(paramName) {
    _getParamAsValue(paramName, false, null)
} 

def optional(paramName) {
    _getParamAsValue(paramName, true, '')
}

def optionalBoolean(paramName) {
   def val = optionalOrDefault(paramName, false)
   return val && !val.toString().trim().equalsIgnoreCase("false")
}

def optionalOrDefault(paramName, defaultValue) {
    _getParamAsValue(paramName, true, defaultValue)
}

def optionalOption(prefix, paramName) {
    def val = _getParamAsValue(prefix + paramName, true, '')
    if (!val) {
        return '';
    } 
    return " --${paramName} ${val} "
}

def optionalOptionWithKey(prefix, optionKey, paramName) {
    def val = _getParamAsValue(prefix + paramName, true, '')
    if (!val) {
        return '';
    } 
    return "--${optionKey} ${val}"
}

def isIonTorrent(reads) {
   try {  
        if (params.containsKey('seq_type')) {
            return (params.seq_type == 'ion')
        }           
        return isCompatibleWithSeqType(reads, 'ion', null)
    } catch(Throwable t) {
        exit 1, "could not get 'seq_type', unexpected exception: ${t.asString()}"
    }     
}

def isNanopore(reads) {
   try {  
        if (params.containsKey('seq_type')) {
            return (params.seq_type == 'nanopore')
        }  
        return isCompatibleWithSeqType(reads, 'nanopore', null)
    } catch(Throwable t) {
        exit 1, "could not get 'seq_type', unexpected exception: ${t.asString()}"
    }     
}

def isIlluminaPaired(reads) {
   try {  
        if (params.containsKey('seq_type')) {
            return (params.seq_type == 'illumina_paired')
        }    
        return isCompatibleWithSeqType(reads, 'illumina_paired', null)    
    } catch(Throwable t) {
        exit 1, "could not get 'seq_type', unexpected exception: ${t.asString()}"
    }     
}

def isCompatibleWithSeqType(reads, types, entrypoint) {
   try {  
        if (!reads || !types) {
            return false
        }
        def allowed = types instanceof java.util.Collection ? types : [ types ]
        def type = params.containsKey('seq_type') ? params.seq_type : ''
        if (!type) {
            if (reads instanceof java.util.Collection && (reads.size() > 1)) {
                type = 'illumina_paired'
            } else {
                def seReads = (reads instanceof java.util.Collection) ? reads[0] : reads
                def gzipped = seReads.toFile().getName() ==~ /.+\.gz/            
                seReads.toRealPath().toFile().withInputStream {
                    try (InputStream is = (gzipped ? new java.util.zip.GZIPInputStream(it) : it);
                        InputStreamReader isr = new InputStreamReader(is);
                        BufferedReader br = new BufferedReader(isr);) {
                            // look at the header to understand if this was generated by nanopore
                            // XXX we need a more robust check
                            /*
                            The headers are defined by the basecaller (guppy, Albacore, bonito, etc) during basecalling.
                            They have changed many times over the years, but generally they start with a readID of @+UUIDv4 
                            followed by multiple entries separated by spaces, of form key=value
                            */
                            type = (br.readLine() ==~ /^@[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}.*/ ) ? 'nanopore' : 'ion'
                    }
                }
            }
        }
        def check = allowed.find { it.trim() == type.trim() }.any()

        if (!check && entrypoint && params.incompatible_step_warning) {
            log.warn "${entrypoint} not compatible with seq type: ${type}"
        }
        return check
    } catch(Throwable t) {
        exit 1, "could not get 'seq_type', unexpected exception: ${t.asString()}"
    }     
}

def isSegmentedMapping() {
   try {  
        // XXX use only parameter 'segmented_mapping' should be used
        return (params.containsKey('segmented_mapping') && params.segmented_mapping) || 
        (workflow.scriptName && workflow.scriptName.endsWith('_segmented.nf'))  
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }     
}

def checkEnum(enumClass, name) {
    if (!enumClass.values().any { it.name() == name }) {
        log.warn("'${name}' is not a valid option. Expecting one of: ${enumClass.values()}")
        return false;
    }
    return true
}

def getHostReference() {
    try {  
        if (!params.containsKey('host')) {
            exit 2, "host reference code should be provided";
        }     
        if (!params.containsKey('hosts_dir')) {
            exit 2, "hosts_dir should be provided";
        }          
        // CROSSING PARAM
        if (!params.containsKey('riscd')) {
            exit 2, "one of: [ds, riscd] param should be provided";
        } 
        def crossValue = parseRISCD(params.riscd).ds         
        def path = "${params.hosts_dir}/${params.host}{.,_}*.f*a"
        return Channel.fromPath(
                path, checkIfExists: true
        )
        .first().map { [ crossValue, params.host, it ] }        
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }     
}

def getLongReads() {
    if (!params.containsKey('long_reads')) {
        exit 2, "long_reads should be provided";
    }   
    def long_reads = params.long_reads
    if (!long_reads.containsKey('cmp') || !long_reads.containsKey('riscd')) {
        exit 2, "long_reads.cmp and long_reads.riscd should be provided";
    }
    return getResult(long_reads.cmp, long_reads.riscd, null)
}

def paramWrap(paramName, wrap) {
    return wrap.replace('{}', "${param(paramName)}")
}

def optWrap(paramName, wrap) {
    def val = optional(paramName)
    if (!val) {
        return ''
    }
    return wrap.replace('{}', "${val}")
}

def getInputFromPaths(key, paths) {
    return Channel.fromPath(
            paths, checkIfExists: params.check_file_existence
    )        
    .collect()
    .map { [ key, it ] }
}
