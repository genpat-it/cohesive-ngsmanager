#!/usr/bin/env python3

import sys


def create_coverage_report(coverage_file, reference, cmp, ds, coverage_import_file):
    cov = open(coverage_file, "r").readlines()

    with open(coverage_import_file, 'w') as cov_import:
        cov_import.write("CMP_ID,SAMPLE_DS,COV,H_COV,NOTE,PERC_IUPAC,PERC_NS,CONSENSUS_LENGTH\n")
        cov_import.write(
            "{},{},{},{},mapping on {},{},{},{}\n".format(cmp, ds, cov[0].rstrip(), cov[1].rstrip(), reference, '','',''))


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.stderr.write(
            "Usage: %s <coverage_file> <reference> <cmp> <ds> "
            "<output_file_import_coverage>" % (
                sys.argv[0]))
        sys.exit(1)
    else:
        create_coverage_report(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
