include { parseMetadata } from './common.nf'

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
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }    
}

def isRunningFromSampleSheet() {
    return false
}

def isVirus(riscd) {
   try {  
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'virus')
        }           
       return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }     
}

def isWNV(riscd) {
   try {  
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'wnv')
        }           
        return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }     
}

def isNGSMG16S(riscd) {
   try {  
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'NGSMG16S')
        }           
        return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }     
}

def isNegativeControl(riscd) {
   try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'NGTVCTRL')
        }    
        return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }     
}

def isNegativeControlSarsCov2(riscd) {
   try {     
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'NGTVCTRLSC2')

        }          
        return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }     
}

def isPositiveControlSarsCov2(riscd) {
   try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'PSTVCTRLSC2')
        }         
        return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }     
}

def isBacterium(riscd) {
 try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'bacterium')
        }
        return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }        
}

def isAmpliseq(riscd) {
 try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'ampliseq')
        }
        return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }        
}

def isSarsCov2(riscd) {
     try {      
        if (params.containsKey('sample_type')) {
            return (params.sample_type == 'sarscov2')
        }         
        return false
    } catch(Throwable t) {
        exit 1, "could not get sample type, unexpected exception: ${t.asString()}"
    }     
}