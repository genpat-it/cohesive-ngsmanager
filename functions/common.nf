import java.time.format.DateTimeFormatter

//@see https://github.com/nextflow-io/nextflow/issues/1694
//https://github.com/nf-core/sarek/blob/a7679b9b5c178351b1e96a3ffe7ee81ddf9aad06/main.nf#L226
EMPTY_FILE = file("${workflow.workDir}/empty_file")
EMPTY_FILE.text = ''

DT_PATTERN = DateTimeFormatter.ofPattern("yyMMdd");
DT_PATTERN_YEAR = DateTimeFormatter.ofPattern("yyyy");


def logHeader(message){
    // Log colors ANSI codes
    def monolog = params.containsKey('monochrome_logs') && params.monochrome_logs
    def c_reset = monolog  ? '' : "\033[0m";
    def c_cyan = monolog ? '' : "\033[0;36m";
    return """    
        ----------------------------------------------------${c_reset}
        ${c_cyan}scriptId${c_reset}=$workflow.scriptId
        ${c_cyan}start${c_reset}=$workflow.start
        ${c_cyan}commandLine${c_reset}=$workflow.commandLine
        ${c_cyan}launchDir${c_reset}=$workflow.launchDir
        ----------------------------------------------------${c_reset}
        ${c_cyan}  ${message} ${c_reset} 
        ----------------------------------------------------${c_reset}
    """.stripIndent()    
}

def flattenPath(filename){
   try {        
       filename.substring(filename.lastIndexOf("/")+1)
   } catch(Throwable t) {
       exit 1, "unexpected exception: ${t.asString()}"
   } 
}

def parseMetadata(input) {
    def filename
    if (input instanceof java.util.Collection) {
        filename = input.flatten()[0].getName()
    } else if (input instanceof java.nio.file.Path) {
        filename = input.toFile().getName()
    } else if (input instanceof java.io.File) {
        filename = input.getName()
    } else {
        filename = input
    }
    return parseMetadataFromFileName(filename)
}

def parseMetadataFromFileName(input){
  try {        
        def filename = (input instanceof java.util.Collection) ? input.flatten()[0].getName() : input

        def md = [:]
        def matcher = (filename =~ /^(DS\d+)\D(DT\d+)_(\d{4})(\.[^\.]+\.\d+\.\d+\.\d+)(_[^_]+)?.*$/)

        if (matcher.matches()) {
            md.ds = matcher.group(1)
            md.dt = matcher.group(2)
            md.anno = matcher.group(3)
            md.cmp = matcher.group(3)+matcher.group(4)
            md.suffix = matcher.group(5) //seems not to work properly
        } else {
            exit 1, "unexpected file name: ${filename}"
        }
        return md
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }     
}

/* execution specific stuff */
def executionMetadata() {
    try {        
        if (params.containsKey('analysis_date')) {
            def matcher = (params.analysis_date =~ /^(\d{2})(\d+)$/)

            if (matcher.matches()) {
                def anno = "20${matcher.group(1)}"
                def dt = "${matcher.group(1)}${matcher.group(2)}"
                return [
                    "anno": anno,
                    "dt": "DT${dt}"
                ]
            } else {
                exit 2, "param 'analysis_date' should match pattern 'yyMMddHHmmSSsss' or 'yyMMdd' (e.g.210524 or 210525093000123)"
            }
        } else {
            return [
                "anno":  workflow.start.format(DT_PATTERN_YEAR),
                "dt": "DT" + workflow.start.format(DT_PATTERN)
            ]
        }
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }       
}

def csv2map(_path, separator){
    try {       
        def path = (_path instanceof java.util.Collection) ? _path.flatten()[0] : _path
        def lines = path.toRealPath().toFile().readLines()
        if (lines.size() != 2) {
            exit 1, "${path} should have exactly 2 lines"
        }

        def headers = lines[0].split(separator)
        def values =  lines[1].split(separator)

        if (headers.length != values.length) {
            exit 1, "headers size not matching values size"
        }

        def md = [:]
        for (int i = 0; i < values.length; i++) {
            md["${headers[i]}"] = "${values[i]}";
        }

        return md;
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }    
}

def extractKey(input) {
    try {        
        def value = (input instanceof java.util.Collection) ? input.flatten()[0] : input

        value = (value instanceof  java.nio.file.Path) ? value.getName() : value

        def matcher = (value =~ /^(DS\d+)(\D.*)?$/)

        if (matcher.matches()) {
            //println "${value} -> ${matcher.group(1)}"
            return matcher.group(1)
        } else {
            matcher = (value =~ /^\d+-(\d+)-[^-]+-.+$/)
            if (matcher.matches()) {
                //println "${value} -> ${matcher.group(1)}"
                return 'DS'+matcher.group(1)
            }
        }
        exit 1, "could not extract DS or RISCD from: ${input}"
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }
}

def extractDsRef(input) {
    try {     
        def value = (input instanceof java.util.Collection) ? input.flatten().find { it instanceof  java.nio.file.Path } : input

        value = (value instanceof  java.nio.file.Path) ? value.getName() : value

        def matcher = (value =~ /^(DS\d+)(\D.*)_([^_]+)\..+$/)

        if (matcher.matches()) {
            // use string, NOT GStringImpl (e.g. "${var}")
            return matcher.group(1) + "-" + matcher.group(3)
        } else {
            exit 1, "could not extract Base Name from: ${input}"
        }
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }
}

def _isRunningFromSampleSheet() {
    return params.containsKey('samplesheet') && (params.samplesheet instanceof CharSequence);
}

def getBaseName(fileName) {
    try {      
        return (fileName - ~/\.\w+$/)
    } catch(Throwable t) {
        exit 1, "could not get getBaseName, unexpected exception: ${t.asString()}"
    }       
}

def getGB(inputPath) {
    try {
        def referencePath = (inputPath instanceof java.util.Collection) ? inputPath.flatten()[0] : inputPath
        def realPath = referencePath.toRealPath()
        def refName = referencePath.getName()
        if (refName ==~ ("^.+\\.gb\$")) {
            // already a GB
            return referencePath
        }
        if (! refName ==~ ("^.+\\.fa[^\\.]*\$")) {
            // unexpected file name
            log.warn "reference supposed to be a FASTA or GB"
            return EMPTY_FILE
        }
        def expectedGBName = refName.replaceFirst("^(.+)\\.fa[^\\.]*\$", "\$1.gb")
        if (refName == expectedGBName) {
            if (refName != EMPTY_FILE.getName()) {
                log.warn "unexpected file name for fasta reference: ${refName}"  
            }
            return EMPTY_FILE
        }
        def result = realPath.resolveSibling(expectedGBName)
        if (!result || !result?.exists()) {
          log.warn "not found GenBank reference to run snippy for: ${refName}"
          return EMPTY_FILE
        }
        return result;
    } catch(Throwable t) {
        log.error "error while looking for a GB reference for: ${inputPath}, exception: ${t.asString()}"
        return EMPTY_FILE
    }        
}

def getEmpty() {
    return EMPTY_FILE;
}

def parseRISCD(riscd){
     try {        
        def md = [:]
        def matcher = (riscd =~ /^(\d+)-(\d+)-(.+)-([^_]+)(_.+)?$/)

        if (matcher.matches()) {
            md.dt = "DT${matcher.group(1)}"
            md.ds = "DS${matcher.group(2)}"
            md.acc = matcher.group(3)
            md.met = matcher.group(4)
            md.ref = matcher.group(5) ? matcher.group(5).substring(1) : ''
        } else {
            exit 1, "unexpected RISCD format: ${riscd}"
        }
        log.debug "RISCD parsing: ${md.toString()}"
        return md
    } catch(Throwable t) {
        exit 1, "unexpected exception: ${t.asString()}"
    }    
}

def taskMemory(obj, attempt) {
    try {
        def mem 
        switch(attempt) {
            case 1: 
                mem = obj; break
            case 2:
                mem = obj * 3; break
            default: 
                mem = params.max_memory as nextflow.util.MemoryUnit
        }
        return (mem.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1) ? params.max_memory : mem
    } catch (t) {
        println "taskMemory: unexpected exception: ${t.asString()}"
        return obj
    }
}
def taskCpus(obj, attempt) {
    try {
        def max_cpus = params.max_cpus as int
        return obj > max_cpus ? max_cpus : obj
    } catch (t) {
        println "taskMemory: unexpected exception: ${t.asString()}"
        return obj
    }
}

def taskTime(obj, attempt) {
    try {
        def time 
        switch(attempt) {
            case 1: 
                time = obj; break
            default: 
                time = params.max_time as nextflow.util.Duration
        }
        return (time.compareTo(params.max_time as nextflow.util.Duration) == 1) ? params.max_time : time
    } catch (t) {
        println "taskTime: unexpected exception: ${t.asString()}"
        return obj
    }
}

def getRisCd(md, ex, step, method) {
  try {
    return "${ex.dt.drop(2)}-${md.ds.drop(2)}-${step}-${method}"
  } catch(Throwable t) {
      exit 1,  "could not get RISCD from input for: ${md}, ${ex}, ${step}, ${method}; error: ${t.message}"
  }
}

def stepInputs(riscd_input, md, ex, step, method, opt) {
    try {
      def result  = [:]
      if (riscd_input instanceof java.util.Collection) {
            riscd_input = riscd_input.findAll { it != '-' && it != null}
            if (riscd_input.size() == 0) {
                log.warn "riscd_input should have at least one value"
                return '{}'
            }
      } 
      result.input = riscd_input
      result.riscd = getRisCd(md, ex, step, method)
      def meta = workflow.toMap().findAll { it.key in ['scriptId', 'scriptName', 'commitId', 'revision', 'profile', 'sessionId'] && it.value != null }  ?: [:]
      ['_meta_commit_id', '_meta_commit_date', '_meta_build_tag', '_meta_pipeline_timestamp'].findAll { params.containsKey(it) }.each { meta[(it -~ /^_meta_/)] = params[it] }
      meta.nextflow = nextflow.toJsonMap().findAll { it.key in ['version', 'build'] }
      result.meta = meta
      result.params = opt ?: [:]
      return groovy.json.JsonOutput.toJson(result)
          } catch(Throwable t) {
        exit 1,  "could not track input for: ${md}, error: ${t}"
    }
}

def isFastqRiscd(riscd) {
    try {
        def acc = parseRISCD(riscd)?.acc;
        return acc ==~ /0SQ.+/ || acc ==~ /1PP.+/
    } catch(Throwable t) {
        exit 1,  "could not execute isFastqRiscd for input: ${riscd}, error: ${t}"
    }
}

def isFastaRiscd(riscd) {
    try {
        def acc = parseRISCD(riscd)?.acc 
        return acc ==~ /2AS.+/
    } catch(Throwable t) {
        exit 1,  "could not execute isFastaRiscd for input: ${riscd}, error: ${t}"
    }
}

def isImportedRiscd(riscd) {
    try {
        def met = parseRISCD(riscd)?.met;
        return met == 'external'
    } catch(Throwable t) {
        exit 1,  "could not execute isFastaRiscd for input: ${riscd}, error: ${t}"
    }
}

def isSpeciesSupported(gsp, whitelist, reads, entrypoint) {
  try {
    def genus_species = gsp ? gsp.toLowerCase() : ''
    def allowedList = whitelist.collect { it.toLowerCase() }
    def (genus, species) = genus_species.contains("_") ? genus_species.split('_') : [ genus_species, null ]
    if (allowedList.contains(genus_species)) {
        // genus_species allowed
        return true;
    }       
    if (allowedList.contains(genus)) {
        // ALL genus allowed
        return true;
    }
    if (entrypoint && reads) {
        def (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, '']
        def md = parseMetadataFromFileName(r1.getName())
        if (workflow.scriptName && workflow.scriptName.startsWith('step_')) {
            log.warn "${entrypoint} [${md?.cmp}] - genus_species not supported: ${genus_species}. Supported: ${allowedList}"
        }
    }  
    return false
  } catch(Throwable t) {
      exit 1, "unexpected exception: ${t.asString()}"
  } 
}

def filename(input) {
    if (input instanceof java.util.Collection) {
        return input.flatten()[0].getName()
    } else if (input instanceof java.nio.file.Path) {
        return input.toFile().getName()
    } else if (input instanceof java.io.File) {
        return input.getName()
    } else {
        // bad move
        return input
    }
}