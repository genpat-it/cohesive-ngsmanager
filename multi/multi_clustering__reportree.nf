nextflow.enable.dsl=2

include { getVCFs;param;optionalOrDefault } from '../functions/parameters.nf'
include { taskMemory;getEmpty } from '../functions/common.nf'

if (getReportreeInputType() == 'alleles') {
  include { multi_clustering__reportree } from "../multi/multi_clustering__reportree_alleles"
  include { getAlleles as inputFn } from "../multi/multi_clustering__reportree_alleles"
} else if (getReportreeInputType() == 'alignment') {
  include { multi_clustering__reportree } from "../multi/multi_clustering__reportree_alignment"
  include { getInput as inputFn } from '../functions/parameters.nf'
} else {
  include { multi_clustering__reportree } from "../multi/multi_clustering__reportree_vcf"
  include { getVCFs as inputFn } from '../functions/parameters.nf'
}

def getReportreeInputType() {
    def res = param('multi_clustering__reportree__input')
    if (!(res in ['alleles', 'vcf', 'alignment'])) {
        exit 2, "params (multi_clustering__reportree__input) not valid"    
    } 
    return res
}

workflow {    
    multi_clustering__reportree(inputFn(),  param('metadata'), param('geodata'), optionalOrDefault('multi_clustering__reportree__nomenclature', getEmpty()));
}