#!/usr/bin/env python3

import os
import sys


def process(ds, blast_result, loci_schema_pasteur, output_file):
    if os.path.isfile(blast_result) and os.path.isfile(loci_schema_pasteur):
        out_temp = "%s_allele_tmp.txt" % ds
        cmd = "grep -v WARNING: %s | grep exact - | cut -d'.' -f2  > %s" % (blast_result, out_temp)
        os.system(cmd)
        with open(output_file, 'w') as res, open(out_temp, 'r') as alleles, open(loci_schema_pasteur, 'r') as loci_schema:
            diz_call = {}
            for a in alleles.readlines():
                lcall, acall = a.strip().split("-")
                diz_call.update({lcall: acall})
            head = "Sample"
            info = ds
            loci_list = loci_schema.readlines()[0].split(",")
            for locus in loci_list:
                call = diz_call.get(locus, "LNF")
                head += "\t" + locus
                info += "\t" + call
            res.write(head + "\n")
            res.write(info + "\n")
        os.remove(out_temp)
    else:
        sys.stderr.write('could not open input files')
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.stderr.write("Usage: python %s <ds> <blast_mlst_result> <loci_schema_pasteur> <output_file>" % (sys.argv[0]))
        sys.exit(1)
    else:
        process(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4])
