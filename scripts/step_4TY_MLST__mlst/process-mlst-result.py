#!/usr/bin/env python3

import string
import sys
import os


def add_cc(input_file, output_file, mlst_db_folder):
    summary = open(input_file, "r").readlines()
    complete = open(output_file, "w")
    complete.write("DS;GENPAT;SCHEMA;CC;ST;loci(varianti)\n")
    for riga in summary:
        sample = riga.split("\t")[0]
        ds = sample.split("DS")[1].split("-")[0]
        genpat = sample.split("/")[-1].split("_")[1]
        schema = riga.split("\t")[1]
        st = riga.split("\t")[2].strip()
        loci = " ".join(riga.split("\t")[3:])
        if schema == '-':
            complete.write("%s;%s;%s;;;%s" % (ds, genpat, schema, loci))
            break;           
        try:
            schema_file = open(mlst_db_folder + "/" + schema + "/" + schema + ".txt", "r").readlines()
            header = schema_file[0]
            cc = ""
            if "Lineage" in header:
                for combination in schema_file[1:]:
                    if combination.split("\t")[0].strip() == st:
                        cc = combination.split("\t")[-2].strip()
            else:
                for combination in schema_file[1:]:
                    if combination.split("\t")[0].strip() == st:
                        cc = combination.split("\t")[-1].strip()
            complete.write("%s;%s;%s;%s;%s;%s" % (ds, genpat, schema, cc, st, loci))
        except Exception as ex:
            print("exception :" + str(ex))
            complete.write("%s;%s;%s;;%s;%s" % (ds, genpat, schema, st, loci))
            with open(os.path.splitext(output_file)[0] + "_exception.log", 'aw') as log:
                log.write(
                    "Error while looking for clonal complex in: %s/%s/%s.txt " % (mlst_db_folder, schema, schema))


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.stderr.write("Usage: %s <input_file> <base_name> <db_folder>" % (sys.argv[0]))
        sys.exit(1)
    else:
        add_cc(sys.argv[1], sys.argv[2], sys.argv[3])
