#!/usr/bin/env python3

import sys
from Bio import SeqIO


def clean_denovo(l200_file, calls_file, reference, output_file):
    reference_plain_name = reference.replace("_", "") 
    abricateout = open(calls_file, "r").readlines()
    node = []
    for riga in abricateout:
        if reference_plain_name in riga.replace("_", ""): #NC_045512 -> NC045512 FIXME it's more robust to modify nextflow steps instead and look for reference
            node.append(riga.split()[1].strip())
    filt_sequences = []
    denovo = open(l200_file, 'r')
    for record in SeqIO.parse(denovo, "fasta"):
        head_seq = record.id
        if head_seq in node:
            filt_sequences.append(record)
    if len(filt_sequences) > 0:
        with open(output_file, "w") as output_handle:
            SeqIO.write(filt_sequences, output_handle, "fasta")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.stderr.write("Usage: %s <input_file> <cmp> <reference> <output_file>" % (sys.argv[0]))
        sys.exit(1)
    else:
        clean_denovo(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
