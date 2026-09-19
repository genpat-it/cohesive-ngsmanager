nextflow.enable.dsl=2

// --- IMPORTS ---
include { step_1PP_trimming__fastp } from '../steps/step_1PP_trimming__fastp'
include { step_3TX_class__kraken2 } from '../steps/step_3TX_class__kraken2'
include { step_2AS_denovo__spades } from '../steps/step_2AS_denovo__spades'
include { getSingleInput } from '../functions/parameters.nf'

// --- GLOBALS ---

// --- INLINE PROCESSES ---

// --- SUB WORKFLOWS ---
workflow READS_PROCESSING {
    take:
        rawreads

    main:
        trimmed_out = step_1PP_trimming__fastp(rawreads)
        trimmed = trimmed_out.trimmed

        assembly_out = step_2AS_denovo__spades(trimmed)
        kraken_out = step_3TX_class__kraken2(trimmed)
    emit:
        assembled = assembly_out.assembled
        genus_report = kraken_out.genus_report

}

// --- ENTRYPOINT ---
workflow {
    rawreads = getSingleInput()
    READS_PROCESSING(rawreads)
}