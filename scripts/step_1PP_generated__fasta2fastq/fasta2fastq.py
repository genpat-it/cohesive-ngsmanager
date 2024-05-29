#!/usr/bin/env python3

import string 
import os
import argparse
from Bio.Seq import Seq

if __name__ == "__main__" :
	parser = argparse.ArgumentParser(description='fasta2fastq')
	parser.add_argument('--fasta', help='input folder', type=str, required=True)
	parser.add_argument('--fastq', type=str, help="output folder", required=True)
	args = parser.parse_args() 
	outdir = args.fastq
	indir = args.fasta
	if os.path.exists(outdir) == False:
		os.system("mkdir "+outdir)

	for fastafile in os.listdir(indir):
		nomefastqR1 = fastafile.replace(".fasta","_R1.fastq")
		R1 = open(outdir+"/"+nomefastqR1,"w")
		nomefastqR2 = fastafile.replace(".fasta","_R2.fastq")
		R2 = open(outdir+"/"+nomefastqR2,"w")
		fasta = open(indir+"/"+fastafile).readlines()
		sequenza = ""
		header = fasta[0][1:].strip()
		for riga in fasta[1:]:
			sequenza = sequenza + riga.strip()
		numeri = range(0,len(sequenza))
		c = 2
		for numero in numeri:
			if numero + 150 < len(sequenza):
				pezzo = sequenza[numero:numero+150]
			else:
				#print str(numero)+":"+str(len(sequenza))
				pezzo = sequenza[numero:len(sequenza)]
			my_dna = Seq(pezzo)
			R1.write("@"+header+"-"+str(c)+"/1\n"+pezzo+"\n+\n"+("I"*len(pezzo))+"\n")
			rc = my_dna.reverse_complement()
			R2.write("@"+header+"-"+str(c)+"/2\n"+str(rc)+"\n+\n"+("I"*len(pezzo))+"\n")
			c+=2
		R1.close()
		R2.close()