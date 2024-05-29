#!/usr/bin/env python3

import os
import argparse

if __name__ == "__main__" :
	parser = argparse.ArgumentParser(description='chewbbaca')
	# parser.add_argument('--outDir',help='output folder',type=str, required=True)
	parser.add_argument('--stat',type=str,help="stats file chewbbaca",required=True)
	args = parser.parse_args() 
	# if args.outDir==".":
	# 	outdir = os.getcwd()
	# else:
	# 	outdir = args.outDir
	if os.path.exists(args.stat):
		s = open(args.stat,'r')
		data = s.readlines()[1]
		header = "genome,calledPerc,calledNum,annotated,new,notFound,discarded"
		g,exc,inf,lnf,plot,niph,alm,asm = data.split("\t")
		total = int(exc)+int(inf)+int(lnf)+int(plot)+int(niph)+int(alm)+int(asm)
		perc = round((float(exc)+float(inf))/total,3)
		called = int(exc)+int(inf)
		disc = int(plot)+int(niph)+int(alm)+int(asm)
		data_out = g+","+str(perc)+","+str(called)+","+exc+","+inf+","+lnf+","+str(disc)
		print(header)
		print(data_out)
	else:
		print("Statistics file not found "+args.stat)
		exit(2)