#!/usr/bin/env python3

from Bio import SeqIO
import argparse


def filter_seq(output_name, check_name, scaffolds, min_length, min_cov):
    filtered_sequences = []
    with (open(scaffolds, 'r') as input_handle,
          open(output_name, "w") as output_handle,
          open(check_name, "w") as check_handle):
        check_handle.write("Nodes_ID\tCHECK\n")

        for record in SeqIO.parse(input_handle, "fasta"):
            head_seq = record.id
            coverage = float(head_seq.split('_')[5])
            if len(record.seq) >= min_length and coverage >= min_cov:
                filtered_sequences.append(record)
                check_handle.write(head_seq + "\tPASS\n")
            else:
                check_handle.write(head_seq + "\tFAILED\n")

        SeqIO.write(filtered_sequences, output_handle, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='FASTQ Report')
    parser.add_argument('-o', '--out', help='Output file name', type=str, required=True)
    parser.add_argument('-oc', '--out_check', help='Output check file name', type=str, required=True)
    parser.add_argument('-f', '--fasta', help='Scaffolds fasta', type=str, required=True)
    parser.add_argument('-l', '--length', help='Min length', type=int, required=True)
    parser.add_argument('-c', '--cov', help='Min kmerCov ', type=int, required=True)
    args = parser.parse_args()

    filter_seq(args.out, args.out_check, args.fasta, args.length, args.cov)
