include { parseMetadata } from './common.nf'

def sortSampleSheet(input, output, batchSize) {
    try {
        def readlines = new File(input).readLines()
        if (readlines.size() > 1) {
            // clone lines (needed later)
            def allLines = readlines.collect()
            
            // assume first row as header
            def headers = readlines[0].split(",")

            // calculate numer of buckets (i.e. batch to be processed)
            def totalLines = readlines.size()
            def bucketsNum = Math.ceil((totalLines -1) / batchSize) as int
            log.info "\ntotalLines (header included): ${totalLines}, batchSize: ${batchSize}, buckets: ${bucketsNum}"
                
            def linesToWrite
            if (bucketsNum > 1) {
                // init buckets
                def buckets = [:]
                (0..<bucketsNum).each{buckets[it] = []}
                
                // remove header
                def header = readlines.remove(0)

                // sort samplesheet by sample type
                def sortedLines = readlines.sort { _getSampleType(headers, it, allLines) }

                // 'smarmella' samples 
                def nextBucket = 0
                (0..<sortedLines.size()).each {
                    buckets[nextBucket++] << sortedLines[it]
                    nextBucket = nextBucket %= bucketsNum;
                }

                // only last batch should have less then 'batchSize' samples
                (0..<(buckets.size()-1)).each {
                    for (int y = buckets[it].size(); y <  batchSize; y++) {
                        buckets[it] << buckets[bucketsNum-1].removeLast()
                    }
                }
                linesToWrite = buckets.values().sum()
                
                // header on top
                linesToWrite.add(0, header)

                log.info "\nBuckets:"
                buckets.values().each { lines -> lines.each { log.info "${it.split(',')[0]} ${it.split(',')[8]}"}; log.info "-" }
            } else {
                linesToWrite = allLines
            }
            assert totalLines == linesToWrite.size()

            new File(output).withWriter {
                linesToWrite.each { line ->
                    it.writeLine(line)
                }
            }
        }
    } catch (t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }
}

def filterSampleSheet() {
    try {
        if (!isRunningFromSampleSheet()) {
            exit 1, "required params (samplesheet) not found";
        }
        if (!params.containsKey('filtered_samplesheet')) {
            exit 1, "required params (filtered_samplesheet) not found";
        }
        def filename = params.samplesheet
        def sampleSheetFile = new File(filename)
        def readlines = sampleSheetFile.readLines()
        def firstSampleIndex = readlines.findIndexOf { it.matches("^\\d.*")}

        def headerIndex = readlines.findIndexOf { it ==~ /^\d.*/ } -1
        if (headerIndex < 0) {
            exit 1, "Cannot find header in samplesheet: filename: ${filename}"
        }
        //guess separator
        def commas = readlines[headerIndex].findAll(",").size()
        def semicolons = readlines[headerIndex].findAll(";").size()
        if (commas < semicolons) {
            //forse usage of ',' as separator
            readlines = readlines.collect {
                it.replaceAll(",", "__").replaceAll(";", ",")
            }
        }
        def headers = readlines[headerIndex].split(",")

        def maxSamples = params.containsKey('max_samples') ? params.max_samples : Integer.MAX_VALUE
        def includes = params.containsKey('includes') ? Arrays.asList(params.includes.split(",")) : []
        def excludes = params.containsKey('excludes') ? Arrays.asList(params.excludes.split(",")) : []
        def validIncludeExcludeValues = ['bacteria', 'viruses', 'sarscov2', 'tender_ecdc','negative_control', 'positive_control_sarscov2', 'negative_control_sarscov2', 'ampliseq']
        if (includes.intersect(validIncludeExcludeValues).size() != includes.size() ||
                excludes.intersect(validIncludeExcludeValues).size() != excludes.size()) {
            exit 1, "wrong includes/excludes param, valid values: ${validIncludeExcludeValues}"
        }

        def filteredLines = []

        def index = 0
        def samplesConsidered = 0
        readlines.each {line ->
            if (samplesConsidered >= maxSamples) {
                return
            }
            if (index++ < firstSampleIndex) {
                filteredLines << line
                return
            }
            if (_isRowToBeConsidered(headers, includes, excludes, line, readlines)) {
                filteredLines << line
                samplesConsidered++
            }
        }
        if (filteredLines.size() > 1) {
            new File(params.filtered_samplesheet).withWriter {
                filteredLines.each { line ->
                    it.writeLine(line)
                }
            }
        }
    } catch (t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }
}

def _getSampleType(headers, line, lines) {
    def dsIndex = headers.findIndexOf { it.equalsIgnoreCase("SAMPLE_ID") }
    if (dsIndex == -1) {
        log.warn "could not find header column SAMPLE_ID"
        return false
    }
    def ds = line.split(",")[dsIndex]
    def metadata = getDsMetadataFromContent(ds, lines)

    if (isMetadataOfTenderECDC(metadata)) {
        return 'tender_ecdc'
    } else if (isMetadataOfNegativeControl(metadata)) {
        return 'negative_control'
    } else if (isMetadataOfPositiveControlSarsCov2(metadata)) {
        return 'positive_control_sarscov2'
    } else if (isMetadataOfNegativeControlSarsCov2(metadata)) {
        return 'negative_control_sarscov2'
    } else if (isMetadataOfSarsCov2(metadata)) {
        return 'sarscov2'
    } else if (isMetadataOfNGSMG16S(metadata)) {
        return 'NGSMG16S'        
    } else if (isMetadataOfWNV(metadata)) {
        return 'wnv'
    } else if (isMetadataOfVirus(metadata)) {
        return 'viruses'
    } else if (isMetadataOfBacterium(metadata)) {
        return 'bacteria'
    } else if (isMetadataOfAmpliseq(metadata)) {
        return 'ampliseq'
    } else {
        return 'generic'
    }
}

def _isRowToBeConsidered(headers, includes, excludes, line, lines) {   
    if (!includes && !excludes) {
        //nothing to filter
        return true
    }
    def dsIndex = headers.findIndexOf { it.equalsIgnoreCase("SAMPLE_ID") }
    if (dsIndex == -1) {
        log.warn "could not find header column SAMPLE_ID"
        return false
    }
    def ds = line.split(",")[dsIndex]
    def metadata = getDsMetadataFromContent(ds, lines)
    def typeList
    if (includes) {
        //process 'includes' first
        typeList = includes
    } else {
        typeList = excludes
    }
    def result = false
    
    def tender = isMetadataOfTenderECDC(metadata)
    def bacterium = isMetadataOfBacterium(metadata)
    def NGSMG16S = isMetadataOfNGSMG16S(metadata)
    def wnv = isMetadataOfWNV(metadata)
    def virus = isMetadataOfVirus(metadata)
    def sarsCov2 = isMetadataOfSarsCov2(metadata)
    def negativeControl = isMetadataOfNegativeControl(metadata)
    def positiveControlSarsCov2 = isMetadataOfPositiveControlSarsCov2(metadata)
    def negativeControlSarsCov2 = isMetadataOfNegativeControlSarsCov2(metadata)
    def ampliseq = isMetadataOfAmpliseq(metadata)

    // have to check multiple conditions sometimes
    if ('tender_ecdc' in typeList) {
        result = result || tender
    }
    if ('sarscov2' in typeList) {
        result = result || (!tender && sarsCov2)
    }
    if ('negative_control' in typeList) {
        result = result || negativeControl
    }
    if ('positive_control_sarscov2' in typeList) {
        result = result || positiveControlSarsCov2
    }
    if ('negative_control_sarscov2' in typeList) {
        result = result || negativeControlSarsCov2
    }
    if ('bacteria' in typeList) {
        result = result || (!tender && !sarsCov2 && !negativeControl && !positiveControlSarsCov2 && !negativeControlSarsCov2 && bacterium)
    }
    if ('viruses' in typeList) {
        result = result || (!tender && !sarsCov2 && !negativeControl && !positiveControlSarsCov2 && !negativeControlSarsCov2 && virus)
    }
    if ('wnv' in typeList) {
        result = result || (!tender && !sarsCov2 && !negativeControl && !positiveControlSarsCov2 && !negativeControlSarsCov2 && wnv)
    }
    if ('NGSMG16S' in typeList) {
        result = result || (!tender && !sarsCov2 && !negativeControl && !positiveControlSarsCov2 && !negativeControlSarsCov2 && NGSMG16S)
    }    
    if ('ampliseq' in typeList) {
        result = result || (!tender && !sarsCov2 && !negativeControl && !positiveControlSarsCov2 && !negativeControlSarsCov2 && ampliseq)
    }    
    return includes ? result : !result
}

def getSampleSheetMetadata(value){
    try {      
        def riscd = (value instanceof java.util.Collection) ? value.flatten()[0] : value
        def matcher = (riscd =~ /\d+-(\d+)-.+$/)
        if (!matcher.matches()) {
            log.warn "cannot extract DS from riscd: ${riscd}"
            return [:];
        }
        return getDsMetadata("DS${matcher.group(1)}")
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }    
}

def getDsMetadata(ds){
    try {

        if (!isRunningFromSampleSheet()) {
            exit 1, "required param (samplesheet) not found";
        }
        def filename = params.samplesheet

        if (!new File(filename).exists()) {
            exit 1, "Could not find SampleSheet data for ds: ${ds}, filename: ${filename}"
        }
        def lines = new File(filename).readLines();

        return getDsMetadataFromContent(ds, lines)
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }
}

def getDsMetadataFromContent(ds, lines){
    try {
        def headerIndex = lines.findIndexOf { it ==~ /^\d.*/ } -1
        if (headerIndex < 0) {
            exit 1, "Cannot find header in samplesheet"
        }

        def headers = lines[headerIndex].split(",")
        if (headers == null) {
            exit 1, "Could not load SampleSheet header data for ds: ${ds}"
        }
        def cmpLine =  lines.find { it.matches(String.format("^[^,]+,%s,.+",ds)) }?.split(",")

        if (!cmpLine) {
            exit 1, "Could not load SampleSheet data for ds: ${ds}"
        }

        if (headers.length < cmpLine.length) {
            log.debug "headers size not matching values size -ds: ${ds}"
        }

        def md = [:]
        for (int i = 0; i < cmpLine.length && i < headers.length; i++) {
            md["${headers[i]}".toString()] = "${cmpLine[i]?.trim()}".toString();
        }

        return md;
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }
}

def isRunningFromSampleSheet() {
    return params.containsKey('samplesheet') && (params.samplesheet instanceof CharSequence);
}

def getRunName() {
     try {        
        if (!isRunningFromSampleSheet())  {
            return "";
        }
        def run_path = params.containsKey('run_path') && params.run_path instanceof CharSequence ? params.run_path : params.samplesheet - ~/\/[^\/]+$/
        def matcher = (run_path =~ /^.+\/([^\/]+)\/?$/)

        if (!matcher.matches()) {
            log.warn "cannot get run name, unexpected run path: ${run_path}"
            return "";
        }
        return matcher.group(1)
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }    
}

def isTenderECDC(riscd) {
   try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'tender_ecdc')
        }       
        if (!isRunningFromSampleSheet())  {
            log.debug "params.samplesheet not set - could not check we are processing a TENDER sample."
            return false
        }
        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfTenderECDC(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isMetadataOfTenderECDC(sampleMetaData) {
   try {      
        def motivo = sampleMetaData.LS_MOTIVO
        if (!motivo) {
            log.warn "column LS_MOTIVO not valorized!"
        }
        return (motivo ==~ /^(.*\D)?201986:.*/)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isVirus(riscd) {
   try {  
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'virus')
        }           
        if (!isRunningFromSampleSheet())  {
            log.debug "params.samplesheet not set - could not check we are processing a VIRUS sample."
            return false
        }

        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfVirus(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isMetadataOfVirus(sampleMetaData) {
   try {      
        def specieCodFamiglia = sampleMetaData.CMP_SPECIE_COD_FAMIGLIA
        if (!specieCodFamiglia) {
            log.warn "column CMP_SPECIE_COD_FAMIGLIA not valorized!"
        }
        def accertCodice = sampleMetaData.ACCERT_CODICE ?: ''
        return (accertCodice in ['G308', ''] && specieCodFamiglia ==~ /^6.*/ && specieCodFamiglia != '60606')        
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isWNV(riscd) {
   try {  
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'wnv')
        }           
        if (!isRunningFromSampleSheet())  {
            log.debug "params.samplesheet not set - could not check we are processing a WNV sample."
            return false
        }

        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfWNV(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isMetadataOfWNV(sampleMetaData) {
   try {      
        def specieCodFamiglia = sampleMetaData.CMP_SPECIE_COD_FAMIGLIA
        if (!specieCodFamiglia) {
            log.warn "column CMP_SPECIE_COD_FAMIGLIA not valorized!"
        }
        def accertCodice = sampleMetaData.ACCERT_CODICE ?: ''
        return (accertCodice in ['G308', 'NGSMGST', ''] && specieCodFamiglia == '60606')        
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isNGSMG16S(riscd) {
   try {  
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'NGSMG16S')
        }           
        if (!isRunningFromSampleSheet())  {
            log.debug "params.samplesheet not set - could not check we are processing a NGSMG16S sample."
            return false
        }

        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfNGSMG16S(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isMetadataOfNGSMG16S(sampleMetaData) {
   try {      
        def accertCodice = sampleMetaData.ACCERT_CODICE ?: ''
        return accertCodice == 'NGSMG16S'     
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isNegativeControl(riscd) {
   try {      
        if (params.containsKey('sample_type')) {
            return false
        }    
        if (!isRunningFromSampleSheet())  {
            log.debug "params.samplesheet not set - could not check we are processing a NEGATIVE CONTROL sample."
            return false
        }
        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfNegativeControl(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isMetadataOfNegativeControl(sampleMetaData) {
   try {      
        def motivo = sampleMetaData.LS_MOTIVO
        if (!motivo) {
            log.warn "column LS_MOTIVO not valorized!"
        }
        def accertCodice = sampleMetaData.ACCERT_CODICE ?: ''
        return (motivo ==~ /^(.*\D)?202135:.*/) || (accertCodice == 'MTGLMBN')
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isNegativeControlSarsCov2(riscd) {
   try {     
        if (params.containsKey('sample_type')) {
            return false
        }          
        if (!isRunningFromSampleSheet())  {
            log.debug "params.samplesheet not set - could not check we are processing a NEGATIVE CONTROL SARSCOV2 sample."
            return false
        }
        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfNegativeControlSarsCov2(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isMetadataOfNegativeControlSarsCov2(sampleMetaData) {
   try {      
        def accertCodice = sampleMetaData.ACCERT_CODICE ?: ''
        return accertCodice in ['MTGSARSN', 'G363']
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}


def isPositiveControlSarsCov2(riscd) {
   try {      
        if (params.containsKey('sample_type')) {
            return false
        }         
        if (!isRunningFromSampleSheet())  {
            log.debug "params.samplesheet not set - could not check we are processing a POSITIVE CONTROL SARSCOV2 sample."
            return false
        }
        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfPositiveControlSarsCov2(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isMetadataOfPositiveControlSarsCov2(sampleMetaData) {
   try {      
        def accertCodice = sampleMetaData.ACCERT_CODICE ?: ''
        return accertCodice in ['MTGSARSP', 'G362']
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isBacterium(riscd) {
 try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'bacterium')
        }
        if (!isRunningFromSampleSheet()) {
            log.debug "params.samplesheet not set - could not check we are processing a BACTERIUM sample."
            return false
        }
        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfBacterium(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }        
}

def isMetadataOfBacterium(sampleMetaData) {
 try {   
        def specieCodFamiglia = sampleMetaData.CMP_SPECIE_COD_FAMIGLIA
        if (!specieCodFamiglia) {
            log.warn "column CMP_SPECIE_COD_FAMIGLIA not valorized!"
        }
        def accertCodice = sampleMetaData.ACCERT_CODICE ?: ''
        return (accertCodice in ['G308', 'G313', ''] && specieCodFamiglia ==~ /^5.*/)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }        
}

def isAmpliseq(riscd) {
 try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'ampliseq')
        }
        if (!isRunningFromSampleSheet()) {
            log.debug "params.samplesheet not set - could not check we are processing an AMPLISEQ sample."
            return false
        }
        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfAmpliseq(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }        
}

def isMetadataOfAmpliseq(sampleMetaData) {
 try {   
        def accertCodice = sampleMetaData.ACCERT_CODICE ?: ''
        return (accertCodice == 'G343')
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }        
}

def isSarsCov2(riscd) {
     try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'sarscov2')
        }         
        if (!isRunningFromSampleSheet()) {
            log.debug "params.samplesheet not set - could not check we are processing a SarsCov2 sample."
            return false
        }
        def sampleMetaData = getSampleSheetMetadata(riscd)
        return isMetadataOfSarsCov2(sampleMetaData)
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def isMetadataOfSarsCov2(sampleMetaData) {
     try {
        if (isMetadataOfPositiveControlSarsCov2(sampleMetaData)) {
            //need this after 'specieCodFamiglia' condition was added
            return false;
        }
        def motivo = sampleMetaData.LS_MOTIVO
        if (!motivo) {
            log.warn "column LS_MOTIVO not valorized!"
        }
        def specieCodFamiglia = sampleMetaData.CMP_SPECIE_COD_FAMIGLIA
        if (!specieCodFamiglia) {
            log.warn "column CMP_SPECIE_COD_FAMIGLIA not valorized!"
        }        
        return motivo ==~ /^(.*\D)?(202020|202118|202119|202127|202143|202163):.*/ || specieCodFamiglia == "61304"
    } catch(Throwable t) {
        exit 1, "could not get sample sheet metadata, unexpected exception: ${t.asString()}"
    }     
}

def getSampleTypeFromDS(riscd) {
    if (isTenderECDC(riscd)) {
        return 'tender_ecdc'
    } else if (isNegativeControl(riscd)) {
        return 'negative_control'
    } else if (isPositiveControlSarsCov2(riscd)) {
        return 'positive_control_sarscov2'
    } else if (isNegativeControlSarsCov2(riscd)) {
        return 'negative_control_sarscov2'
    } else if (isSarsCov2(riscd)) {
        return 'sarscov2'
    } else if (isNGSMG16S(riscd)) {
        return 'NGSMG16S'        
    } else if (isWNV(riscd)) {
        return 'wnv'
    } else if (isVirus(riscd)) {
        return 'viruses' 
    } else if (isBacterium(riscd)) {
        return 'bacteria'
    } else if (isAmpliseq(riscd)) {
        return 'ampliseq'
    } else {
        return 'generic'
    }    
}

def logSampleMetadata(reads, category) {
    try {    
        def prefix = params.containsKey('report_prefix') ? "${params.report_prefix}" : ''
        def fileName = "${prefix}summary_${category}.csv"
        reads.collectFile(name: fileName, storeDir: params.tracedir, newLine : true) {
            def md = parseMetadata(it[1])
            def type = getSampleTypeFromDS(it[0])
            def r1 = (it[1] instanceof java.util.Collection) ? it[1][0] : it[1]
            def count = r1.countFastq()
            if (count < (params.raw_reads_threshold ?: 0)) {
                log.warn "[${type}][${md?.cmp}][${md?.ds}][${category}][${count}]"
            } else {
                log.info "[${type}][${md?.cmp}][${md?.ds}][${category}][${count}]"
            }
            return "${type},${md?.cmp},${md?.ds},${category},${count}"   
        }
    } catch(Throwable t) {
        log.warn "could not log sample metadata, unexpected exception: ${t.asString()}"
    }    
}