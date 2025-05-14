nextflow.enable.dsl=2

include { getMockRawReadsFromSampleSheet;getInput } from '../functions/parameters.nf'
include { getEmpty;extractKey;parseMetadataFromFileName } from '../functions/common.nf'
include { isNegativeControl;isNegativeControlSarsCov2;isPositiveControlSarsCov2;isVirus;isBacterium;isSarsCov2;isTenderECDC;isAmpliseq;isRunningFromSampleSheet; getDsMetadata } from '../functions/samplesheet.nf'

workflow pipeline_empty {
    main:
      if (isRunningFromSampleSheet()) {
        rawReads =  getMockRawReadsFromSampleSheet()  
      } else {
        rawReads =  getInput()
      }

      rawReads.branch {
          tender_ecdc: isTenderECDC(it)
          other: true
      }
      .set { rawReadsBranched }

      rawReadsBranched.other.branch {
          bacteria: isBacterium(it)
          sarscov2: isSarsCov2(it) || isNegativeControlSarsCov2(it) || isPositiveControlSarsCov2(it)
          negative_control: isNegativeControl(it)
          viruses: isVirus(it)
          ampliseq: isAmpliseq(it)
          other: true
      }
      .set { branched }

       /* logs how the samples would be processed */
        rawReadsBranched.tender_ecdc
          .map { 
          def md = parseMetadataFromFileName(it[1])
          log.info "[TENDER-ECDC][${md?.cmp}][${md?.ds}] TENDER ECDC EFSA BRANCH"
        }
        branched.bacteria
          .map {
          def md = parseMetadataFromFileName(it[1])
          log.info "[BACTERIA][${md?.cmp}][${md?.ds}]  BACTERIA BRANCH"
        }
        branched.sarscov2
          .map { 
          def md = parseMetadataFromFileName(it[1])
          log.info "[SARSCOV2][${md?.cmp}][${md?.ds}]  SARSCOV2 BRANCH"
        }
        branched.viruses
          .map { 
          def md = parseMetadataFromFileName(it[1])
          log.info "[VIRUSES][${md?.cmp}][${md?.ds}]  VIRUSES BRANCH"
        }
        branched.negative_control
          .map { 
          def md = parseMetadataFromFileName(it[1])
          log.info "[NEGATIVE-CONTROL][${md?.cmp}][${md?.ds}]  NEGATIVE CONTROL BRANCH"
        }
        branched.ampliseq
          .map { 
          def md = parseMetadataFromFileName(it[1])
          log.info "[AMPLISEQ][${md?.cmp}][${md?.ds}]  AMPLISEQ BRANCH"
        }        
        branched.other.map { 
          def md = parseMetadataFromFileName(it[1])
          log.info "[GENERIC][${md?.cmp}][${md?.ds}] GENERIC BRANCH"
        }
}

workflow {
    pipeline_empty()
}