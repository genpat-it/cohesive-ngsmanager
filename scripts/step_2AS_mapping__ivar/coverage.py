#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC


def create_coverage_report(coverage_file, consensus_file, reference, cmp, ds, output_file, coverage_import_file):
    cov = open(coverage_file, "r").readlines()

    with open(output_file, 'w') as check, open(coverage_import_file, 'w') as cov_import:
        check.write("#Mapping\nReference\tCov\tHCov\n")
        check.write(reference + "\t" + str(cov[0].rstrip()) + "\t" + str(cov[1].rstrip()) + "\n")

        n_number = 0
        iupac_number = 0
        seq_length = 0

        consensus = SeqIO.index(consensus_file, "fasta")
        for seqID in consensus.keys():
            seq = consensus[seqID].seq.upper()
            seq_length = len(seq)
            check.write("\n#Consensus\nTotal_Length: " + str(seq_length) + "\n")
            check.write("\nIUPAC\tCOUNT\tPERC\n")
            for letter in IUPAC.ambiguous_dna.letters:
                if letter not in IUPAC.unambiguous_dna.letters:
                    count = seq.count(letter)
                    if seq_length != 0:
                        perc = round(float(seq.count(letter)) / float(seq_length) * 100, 3)
                    else:
                        perc = 0
                    check.write(letter + "\t" + str(count) + "\t" + str(perc) + "\n")
                    if letter == 'N':
                        n_number = count
                    else:
                        iupac_number += count
            break  # considering only first key!

        perc_ns = round(n_number / float(seq_length) * 100, 3) if seq_length > 0 else 0
        perc_iupac = round(iupac_number / float(seq_length) * 100, 3) if seq_length > 0 else 0

        cov_import.write("CMP_ID,SAMPLE_DS,COV,H_COV,NOTE,PERC_IUPAC,PERC_NS,CONSENSUS_LENGTH\n")
        cov_import.write(
            "{},{},{},{},mapping on {},{},{},{}\n".format(cmp, ds, cov[0].rstrip(), cov[1].rstrip(), reference, perc_iupac,
                                                       perc_ns,seq_length))


if __name__ == "__main__":
    if len(sys.argv) < 8:
        sys.stderr.write(
            "Usage: %s <coverage_file> <consensus_file> <reference> <cmp> <ds> <output_file> "
            "<output_file_import_coverage>" % (
                sys.argv[0]))
        sys.exit(1)
    else:
        create_coverage_report(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],
                               sys.argv[7])
