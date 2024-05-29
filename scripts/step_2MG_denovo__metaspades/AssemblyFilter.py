#!/usr/bin/env python3

from Bio import SeqIO
from numpy import median,average
import sys
import argparse

 
def filt_seq(sampleName, scaffolds, minLength, minCov):
	
	filt_sequences = []
	input_handle=open(scaffolds,'r')
	output_handle = open(sampleName+"_spades_scaffolds_L200.fasta", "w")
	check_handle = open(sampleName+"_spades.check", "w")
	check_handle.write("Nodes_ID\tCHECK\n")

	for record in SeqIO.parse(input_handle, "fasta"):
		headSeq = record.id
		coverage = float(headSeq.split('cov_')[1])
		if len(record.seq) >= minLength and coverage >= minCov:
			filt_sequences.append(record)
			check_handle.write(headSeq+"\tPASS\n")
		else:
			check_handle.write(headSeq+"\tFAILED\n")
	
	SeqIO.write(filt_sequences, output_handle, "fasta")
	input_handle.close()
	output_handle.close()
	check_handle.close()
	return

def filt_seq_unicycler(sampleName, scaffolds, minLength, minCov):
	
	filt_sequences = []
	input_handle=open(scaffolds,'r')
	output_handle = open(sampleName+"_unicycler_scaffolds_L200.fasta", "w")
	check_handle = open(sampleName+"_unicycler.check", "w")
	check_handle.write("Nodes_ID\tCHECK\n")

	for record in SeqIO.parse(input_handle, "fasta"):
		headSeq = record.description
		print(headSeq)
		coverage = float(headSeq.split(' ')[2].lstrip("depth=").strip("x"))
		if len(record.seq) >= minLength and coverage >= minCov:
			filt_sequences.append(record)
			check_handle.write(headSeq+"\tPASS\n")
		else:
			check_handle.write(headSeq+"\tFAILED\n")
	
	SeqIO.write(filt_sequences, output_handle, "fasta")
	input_handle.close()
	output_handle.close()
	check_handle.close()
	return

if __name__=="__main__":
	
	parser = argparse.ArgumentParser(description='FASTQ Report')
	parser.add_argument('-n', '--name',help='Sample name',type=str, required=True)
	parser.add_argument('-f', '--fasta',help='Scaffolds fasta',type=str, required=True)
	parser.add_argument('-l', '--length',help='Min length',type=int, required=True)
	parser.add_argument('-c', '--cov',help='Min kmerCov ',type=int, required=True)
	parser.add_argument('-u', '--unicycler',help='Unicycler fasta',action='store_true')
	args = parser.parse_args()
	
	if not args.unicycler:
		filt_seq(args.name,args.fasta,args.length,args.cov)
	else:
		filt_seq_unicycler(args.name,args.fasta,args.length,args.cov)