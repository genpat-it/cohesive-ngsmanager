nextflow.enable.dsl=2

include { module_reads_processing } from '../modules/module_reads_processing'
include { module_wgs_bacteria } from '../modules/module_wgs_bacteria'
include { module_typing_bacteria } from '../modules/module_typing_bacteria'
include { module_covid_emergency } from '../modules/module_covid_emergency'
include { module_vdraft } from '../modules/module_vdraft'
include { module_westnile } from '../modules/module_westnile'
include { module_qc_quast } from '../modules/module_qc_quast'
include { step_4TY_lineage__pangolin } from '../steps/step_4TY_lineage__pangolin'


include { getSingleInput;param } from '../functions/parameters.nf'
include { getEmpty;extractKey;isFastqRiscd;isFastaRiscd } from '../functions/common.nf'
include { isNegativeControl;isNegativeControlSarsCov2;isPositiveControlSarsCov2;isVirus;isWNV;isBacterium;isSarsCov2;isAmpliseq } from '../functions/sampletypes.nf'

workflow ngsmanager_fastq {
    take: 
      rawReads  
    main:
      try {       
        trimmed = module_reads_processing(rawReads).trimmed_with_data

        trimmed.branch {
            bacteria: isBacterium(it)
            sarscov2: isSarsCov2(it) || isNegativeControlSarsCov2(it) || isPositiveControlSarsCov2(it)
            negative_control: isNegativeControl(it)
            viruses: isVirus(it)
            wnv: isWNV(it)
            ampliseq: isAmpliseq(it)
            other: true
        }
        .set { branched }

        /* bacteria stuff here */
        assemblyBacteria = module_wgs_bacteria(branched.bacteria)
        crossedBacteriaData = branched.bacteria.cross(assemblyBacteria) { extractKey(it) }.multiMap { 
          trimmed: it[0]
          assembly: it[1]
        }
        module_typing_bacteria(crossedBacteriaData.trimmed, crossedBacteriaData.assembly)

        /* sarscov2 stuff here */
        module_covid_emergency(branched.sarscov2)

        /* viruses stuff here */
        // no extra modules executed
        
        /* wnv stuff here */
        module_westnile(branched.wnv)
      } catch (t) {
          exit 1, "unexpected exception: ${t.asString()}"
      }
     
}

workflow ngsmanager_fasta {
    take: 
      fasta
    main:
      try {      
        module_qc_quast(fasta)

        fasta.branch {
            bacteria: isBacterium(it)
            sarscov2: isSarsCov2(it) || isNegativeControlSarsCov2(it) || isPositiveControlSarsCov2(it)
            negative_control: isNegativeControl(it)
            viruses: isVirus(it)
            wnv: isWNV(it)
            ampliseq: isAmpliseq(it)
            other: true
        }
        .set { branched }

        /* bacteria stuff here */
        module_typing_bacteria(Channel.empty(), branched.bacteria)

        /* sarscov2 stuff here */
        step_4TY_lineage__pangolin(branched.sarscov2)       

        /* wnv stuff here */
        module_westnile(branched.wnv)        
      } catch (t) {
          exit 1, "unexpected exception: ${t.asString()}"
      }     
}

workflow pipeline_ngsmanager {
  if (isFastqRiscd(param('riscd'))){
    ngsmanager_fastq(getSingleInput())
  } else if (isFastaRiscd(param('riscd'))){
    ngsmanager_fasta(getSingleInput())
  } else {
    exit 2, "unexpected riscd provided: ${param('riscd')}"
  }
}

workflow {
  pipeline_ngsmanager()
}