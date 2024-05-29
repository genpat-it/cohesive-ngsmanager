#!/usr/bin/env python3

import sys
import argparse
import os

def tsv2diz(fileStat):
	diz = {}
	with open(fileStat,'r') as checkIn:
		for line in checkIn:
			key, val = line.strip().split("\t")
			diz[key]= val
	return diz

def makeReport(name,stat):
	with open(name+"_readsCheck.csv",'w') as checkOut:
		noVal = ",".join(["NULL"]*11)
		checkOut.write("#Sample_name,Total_rawReads,Total_rawMbases,Mean_rawLength,Q30_rawReads,Mean_rawAvgQual,Total_trimReads,trimDiscarded,trimUnpaired,trimPaired,trimPairedMbases,Mean_trimPairedLength,Q30_trimPairedReads,Mean_trimPairedAvgQual,Overlap_trimPaired,WARNING,FAIL\n")
		checkOut.write(name+","+str(stat["number_of_reads"])+","+str(stat["number_of_bases"])+","+str(stat["mean_read_length"])+",NULL,"+str(stat["mean_qual"])+","+noVal+"\n")
	return

# MAIN
if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='FASTQ Report')
	parser.add_argument('-n', '--name',help='Sample name',type=str, required=True)
	parser.add_argument('-s', '--stats',help='Output Nanoplot (NanoStas.txt)',type=str, required=True)
	args = parser.parse_args()

	if os.path.exists(args.stats):
		statDiz = tsv2diz(args.stats)
		makeReport(args.name,statDiz)
	else:
		print("Error, statistic file not found")
	
	
	
	
