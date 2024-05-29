#!/usr/bin/env python2

import os
import sys
import string
import csv
from operator import *

# INFO JSON
dbBacteria=sys.argv[2]+"/Bacteria/DB_summary_complete_genome.txt"
dbViral=sys.argv[2]+"/Viral/DB_summary_complete_genome.txt"

# INPUT
infilenameB = sys.argv[1]
basename = infilenameB.replace("_kmerfinder_bacterial.tsv", "")
outname = basename+"_kmerfinder.check"
infilenameV = basename+"_kmerfinder_viral.tsv"

# Creo diz per parsare le info dbBacteria
dizBacteria={}
for line in open(dbBacteria,'r'):
	sline = line.strip().split("\t")
	key = sline[0]
	dizBacteria.update({key:sline})

# Creo diz per parsare le info dbViral
dizViral={}
for line in open(dbViral,'r'):
	sline = line.strip().split("\t")
	key = sline[0]
	dizViral.update({key:sline})

guessFileBacteria = open(infilenameB, 'r')
guessFileViral = open(infilenameV, 'r')
outfile = open(outname, 'w')

linesBacteria = guessFileBacteria.readlines()
linesViral = guessFileViral.readlines()

# Info Templato Virus 
try:
	genomeViral = linesViral[1].split("\t")[0].strip()
	qcovViral = float(linesViral[1].split("\t")[5].strip())
	tcovViral = float(linesViral[1].split("\t")[6].strip())
	DBinfoViral = dizViral.get(genomeViral)
	assAccViral=DBinfoViral[0]
	guessViral=DBinfoViral[1]
except:
	guessViral = "NoVirus"
	genomeViral="NA"
	qcovViral = "NA"
	tcovViral = "NA"
	assAccViral="NA"

# Info Templato bacteria
try:
	genomeBacteria = linesBacteria[1].split("\t")[0].strip()
	qcovBacteria = float(linesBacteria[1].split("\t")[5].strip())
	tcovBacteria = float(linesBacteria[1].split("\t")[6].strip())
	DBinfoBacteria = dizBacteria.get(genomeBacteria)
	guessBacteria=DBinfoBacteria[2]
	assAccBacteria=DBinfoBacteria[0]
	desBacteria=DBinfoBacteria[1]
except:
	guessBacteria = "NoBacteria"
	genomeBacteria = "NA"
	qcovBacteria = "NA"
	tcovBacteria = "NA"
	assAccBacteria="NA"
	desBacteria= "NA"

# Verifico copertura templato Viral
checkViral=""
if guessViral != "NoVirus":	
	if qcovViral > 50 and tcovViral > 50:
		checkViral="PASS"
	elif qcovViral > 40 or tcovViral > 40:
		checkViral="WARNING"
	else:
		checkViral="ALERT"
else:
	checkViral="FAIL"

# Verifico copertura templato Bacteria
checkBacteria=""
if guessBacteria != "NoBacteria":	
	if qcovBacteria > 50 and tcovBacteria > 50:
		checkBacteria="PASS"
	elif qcovBacteria > 40 or tcovBacteria > 40:
		checkBacteria="WARNING"
	else:
		checkBacteria="ALERT"
else:
	checkBacteria="FAIL"

#print(guessViral+"="+checkViral+":"+guessBacteria+"="+checkBacteria)

# Set species/folder name
species=""
if (checkBacteria in ["PASS","WARNING"]) or (checkViral in ["PASS","WARNING"]):
	if (checkBacteria in ["FAIL","ALERT"]) and (checkViral in ["PASS","WARNING"]):
		species=guessViral
	elif (checkBacteria in ["PASS","WARNING"]) and (checkViral in ["FAIL","ALERT"]):
		species=guessBacteria
	elif checkBacteria=="PASS" and checkViral=="WARNING":
		species=guessBacteria
	elif checkBacteria=="WARNING" and checkViral=="PASS":
		species=guessViral
	elif checkBacteria==checkViral=="WARNING" or checkBacteria==checkViral=="PASS":
		species="Contamination"
else:
	species="NoSpecies"

outfile.write("# filename\tspeciesAssigned\tviralGuess\tcheckViral\tsampleKmerCovViral\ttemplateKmerCovViral\tassembly_accViral\tbacteriaGuess\tcheckBacteria\tsampleKmerCovBacteria\ttemplateKmerCovBacteria\tassembly_accBacteria\tdescBacteria\n")
outfile.write(basename+"\t"+species+"\t"+guessViral+"\t"+checkViral+"\t"+str(qcovViral)+"\t"+str(tcovViral)+"\t"+assAccViral+"\t"+guessBacteria+"\t"+checkBacteria+"\t"+str(qcovBacteria)+"\t"+str(tcovBacteria)+"\t"+assAccBacteria+"\t"+desBacteria+"\n")
outfile.close()

#print(species)

if os.path.exists(species) == False:
	os.mkdir(species)

os.system("mv %s* %s" % (basename, species))
os.system("cp %s/%s ." % (species, outname))	




