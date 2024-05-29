#!/usr/bin/env python3

import argparse
import os


def get_info(data, num):
    try:
        d = data[num]
        d_info = d.strip().split("\t")
        return d_info[0] + "," + d_info[6] + "," + d_info[5]
    except:
        return "null,null,null"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create taxa file for import')
    parser.add_argument('-ds', '--ds', help='DS sample', type=str, required=True)
    parser.add_argument('-o', '--outfile', help='Outfile name', type=str, required=True)
    parser.add_argument('-g', '--genus', help='Braken genus', type=str, required=True)
    parser.add_argument('-sp', '--species', help='Braken species', type=str, required=True)
    parser.add_argument('-k', '--kraken', help='Kraken report', type=str, required=True)
    args = parser.parse_args()

    if os.path.exists(args.genus) and os.path.exists(args.species) and os.path.exists(args.kraken):
        info = []
        with open(args.outfile, 'w') as out:
            out.write(
                "SAMPLE,UNCLASS_PC_READS,UNCLASS_NUM_READS,FIRSTSPE_TAXA,FIRSTSPE_PC_READS,FIRSTSPE_NUM_READS,"
                "FIRSTSPE_COVERAGE,GEN1RANK_TAXA,GEN1RANK_PC_READS,GEN1RANK_NUM_READS,GEN2RANK_TAXA,"
                "GEN2RANK_PC_READS,GEN2RANK_NUM_READS,GEN3RANK_TAXA,GEN3RANK_PC_READS,GEN3RANK_NUM_READS,"
                "SPE2RANK_TAXA,SPE2RANK_PC_READS,SPE2RANK_NUM_READS,SPE3RANK_TAXA,SPE3RANK_PC_READS,"
                "SPE3RANK_NUM_READS,RAW_COVERAGE,TRIMMED_COVERAGE\n") 
            info.append(args.ds.replace("DS", ""))
            # UNCL
            uncl = open(args.kraken).readlines()[0]
            uncl_info = uncl.split("\t")
            info.append(uncl_info[0])
            info.append(uncl_info[1])
            # SPECIES + GENUS
            genus = open(args.genus).readlines()
            species = open(args.species).readlines()
            info.append(get_info(species, 1))
            # Set species coverage to null
            info.append("null")
            info.append(get_info(genus, 1))
            info.append(get_info(genus, 2))
            info.append(get_info(genus, 3))
            info.append(get_info(species, 2))
            info.append(get_info(species, 3))
            # Set coverage to null
            info.append("null,null\n")
            out.write(",".join(info))
    else:
        print("ERROR: input file not exists")
