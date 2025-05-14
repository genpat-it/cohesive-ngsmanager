#!/usr/bin/env python3

from Bio import SeqIO
from numpy import median,average
import sys
import argparse
#import multiprocessing as mp

def safe_open(file, mode='rt'):
	# safe open file
	if file.endswith('.gz'):
		import gzip
		return gzip.open(file, mode)
	else:
		return open(file, mode)
	
def readsInfo(fastq):
	res = {}
	lengthReads=[]
	total_bases = 0
	q30_bases = 0
	avgQual=[]
	handle = safe_open(fastq)
	for record in SeqIO.parse(handle, "fastq"):
		lengthReads.append(len(record.seq))
		phredQualityList = record.letter_annotations['phred_quality']
		avgQual.append(average(phredQualityList))
		total_bases += len(phredQualityList)
		q30_bases += sum(q >= 30 for q in phredQualityList)
	handle.close()
	
	q30Reads = round(q30_bases / float(total_bases) * 100, 2)
	Mbases = round(float(sum(lengthReads))/1000000,2)

	res.update({"lengthReads":lengthReads})
	res.update({"avgQual":avgQual})
	res.update({"q30Reads":q30Reads})
	res.update({"Mbases":Mbases})
	return res

# Class Definition
class SampleRaw:
	
	def __init__(self, read1):
		print("Processing raw reads...")
		self.R1 = readsInfo(read1)
		self.warnings = []
		self.fails = []
	
	def get_R1(self):
		return self.R1
		
	def get_warnings(self):
		return self.warnings
		
	def get_fails(self):
		return self.fails
	
	def getNR(self):
		lengthReads = self.R1.get("lengthReads")
		return len(lengthReads)
	
	def get_q30Reads(self):
		readsQual30 = self.R1.get("q30Reads")
		return readsQual30
		
	def get_Qavg(self):
		avgQual = self.R1.get("avgQual")
		return round(median(avgQual),2)
	
	def get_stats(self):
		lengthReads = self.R1.get("lengthReads")
		avgQual = self.R1.get("avgQual")
		Mbases = round(float(sum(lengthReads))/1000000,2)
		return [len(lengthReads),Mbases,min(lengthReads),max(lengthReads),round(average(lengthReads),2),self.get_q30Reads(),min(avgQual),max(avgQual),round(median(avgQual),2)]
	
	def checkQ30(self):
		q30R = self.get_q30Reads()
		if float(q30R) >= 25:
			return 
		elif float(q30R) >= 5 and float(q30R) < 25:
			self.warnings.append("WARN_rawQ30")
		else:
			self.fails.append("FAIL_rawQ30")
			
	def checkQavg(self):
		avgQualR = self.get_Qavg()
		if float(avgQualR) >= 22:
			return 
		elif float(avgQualR) >= 16 and float(avgQualR) < 22:
			self.warnings.append("WARN_rawQavg")
		else:
			self.fails.append("FAIL_rawQavg")
	
	def checkNR_Bact(self):
		readsTot = self.getNR()
		if float(readsTot) >= 350000:
			return 
		elif float(readsTot) >= 250000 and float(readsTot) < 350000:
			self.warnings.append("WARN_rawNR_bact")
		else:
			self.fails.append("FAIL_rawNR_bact")
			
	def checkNR_Vir(self):
		readsTot = self.getNR()
		if float(readsTot) >= 100000:
			return 
		elif float(readsTot) >= 500000 and float(readsTot) < 100000:
			self.warnings.append("WARN_rawNR_vir")
		else:
			self.fails.append("FAIL_rawNR_vir")

class SampleTrimmed:
	
	def __init__(self, read1):
		print("Processing trimmed reads...")
		self.R1 = readsInfo(read1)
		self.warnings = []
		self.fails = []
	
	def get_R1(self):
		return self.R1
	
	def get_warnings(self):
		return self.warnings
		
	def get_fails(self):
		return self.fails
	
	def getNR_total(self):
		lengthReads = self.R1.get("lengthReads")
		return len(lengthReads)
	
	def getNR(self):
		lengthReads = self.R1.get("lengthReads")
		return len(lengthReads)
	
	def get_q30Reads(self):
		readsQual30 = self.R1.get("q30Reads")
		return readsQual30
		
	def get_Qavg(self):
		avgQual = self.R1.get("avgQual")
		return round(median(avgQual),2)
	
	def get_stats(self):
		lengthReads = self.R1.get("lengthReads")
		avgQual = self.R1.get("avgQual")
		Mbases = round(float(sum(lengthReads))/1000000,2)
		return [self.getNR_total(),0,0,Mbases,min(lengthReads),max(lengthReads),round(average(lengthReads),2),self.get_q30Reads(),min(avgQual),max(avgQual),round(median(avgQual),2)]
	
	def checkQ30(self):
		q30R = self.get_q30Reads()
		if float(q30R) >= 35:
			return 
		elif float(q30R) >= 10 and float(q30R) < 35:
			self.warnings.append("WARN_trimQ30")
		else:
			self.fails.append("FAIL_trimQ30")
			
	def checkQavg(self):
		avgQualR = self.get_Qavg()
		if float(avgQualR) >= 24:
			return 
		elif float(avgQualR) >= 18 and float(avgQualR) < 24:
			self.warnings.append("WARN_trimQavg")
		else:
			self.fails.append("FAIL_trimQavg")
	
	def checkNR_Bact(self):
		readsTot = self.getNR()
		if float(readsTot) >= 300000:
			return 
		elif float(readsTot) >= 50000 and float(readsTot) < 300000:
			self.warnings.append("WARN_trimNR_bact")
		else:
			self.fails.append("FAIL_trimNR_bact")
			
	def checkNR_Vir(self):
		readsTot = self.getNR()
		if float(readsTot) >= 90000:
			return 
		elif float(readsTot) >= 40000 and float(readsTot) < 90000:
			self.warnings.append("WARN_trimNR_vir")
		else:
			self.fails.append("FAIL_trimNR_vir")
	
class Sample:
	
	def __init__(self, sampleName, R1, T1=None):
		self.name = sampleName
		self.Raw = SampleRaw(R1)
		if T1 != None:
			self.Trimmed = SampleTrimmed(T1)
		else:
			self.Trimmed = None
		self.warnings = []
		self.fails = []
	
	def get_name(self):
		return self.name

	def raw_check(self):
		self.Raw.checkQ30()
		self.Raw.checkQavg()
		self.Raw.checkNR_Bact()
		self.Raw.checkNR_Vir()
		with open(self.get_name()+"_SRC_raw.csv",'w') as checkOut:
			results = ",".join(str(x) for x in self.Raw.get_stats())
			checkOut.write("#Sample_name,Total_reads,Total_Mbases,Min_length,Max_length,Mean_length,Q30_reads,Min_avgQual,Max_avgQual,Med_avgQual\n")
			checkOut.write(self.get_name()+","+results+"\n")
		return
	
	def raw_esito(self):
		resWarn = ":".join(self.Raw.get_warnings())
		resFail = ":".join(self.Raw.get_fails())
		with open(self.get_name()+"_SRC_raw.check",'w') as resOut:
			resOut.write("#Sample_name,WARNING,FAIL\n")
			resOut.write(self.get_name()+","+resWarn+","+resFail+"\n")
		return

	def trimm_check(self):
		self.Trimmed.checkQ30()
		self.Trimmed.checkQavg()
		self.Trimmed.checkNR_Bact()
		self.Trimmed.checkNR_Vir()
		with open(self.get_name()+"_SRC_treads.csv",'w') as checkOut:
			results = ",".join(str(x) for x in self.Trimmed.get_stats())
			checkOut.write("#Sample_name,Total_reads,Unpaired_reads,Paired_reads,Paired_Mbases,Min_PairedLength,Max_PairedLength,Mean_PairedLength,Q30_PairedReads,Min_PairedAvgQual,Max_PairedAvgQual,Mean_PairedAvgQual\n")
			checkOut.write(self.get_name()+","+results+"\n")
		return
			
	def trimm_esito(self):
		resWarn = ":".join(self.Trimmed.get_warnings())
		resFail = ":".join(self.Trimmed.get_fails())
		with open(self.get_name()+"_SRC_treads.check",'w') as resOut:
			resOut.write("#Sample_name,WARNING,FAIL\n")
			resOut.write(self.get_name()+","+resWarn+","+resFail+"\n")
		return
	
	def trimDiscarded(self):
		discarded = self.Raw.getNR()-self.Trimmed.getNR_total()
		return discarded
		
	def trimOverlap(self):
		overlap = round(float(self.Pear.getNR_Merged())/float(self.Trimmed.getNR())*100,2)
		return overlap
		
			
	def makeReport(self):
		with open(self.get_name()+"_readsCheck.csv",'w') as checkOut:
			resRaw = self.Raw.get_stats()
			if self.Trimmed != None:
				resTrim = self.Trimmed.get_stats()
				resWarn = ":".join(self.Raw.get_warnings()+self.Trimmed.get_warnings())
				resFail = ":".join(self.Raw.get_fails()+self.Trimmed.get_fails())
				results = ",".join(str(x) for x in [resRaw[0],resRaw[1],resRaw[4],resRaw[5],resRaw[8],resTrim[0],self.trimDiscarded(),resTrim[1],resTrim[2],resTrim[3],resTrim[6],resTrim[7],resTrim[10],0])
			else:
				resTrim = [0,0,0,0,0,0,0,0,0,0,0]
				resWarn = ":".join(self.Raw.get_warnings())
				resFail = ":".join(self.Raw.get_fails())
				results = ",".join(str(x) for x in [resRaw[0],resRaw[1],resRaw[4],resRaw[5],resRaw[8],resTrim[0],0,resTrim[1],resTrim[2],resTrim[3],resTrim[6],resTrim[7],resTrim[10],0])
			checkOut.write("#Sample_name,Total_rawReads,Total_rawMbases,Mean_rawLength,Q30_rawReads,Mean_rawAvgQual,Total_trimReads,trimDiscarded,trimUnpaired,trimPaired,trimPairedMbases,Mean_trimPairedLength,Q30_trimPairedReads,Mean_trimPairedAvgQual,Overlap_trimPaired,WARNING,FAIL\n")
			checkOut.write(self.get_name()+","+results+","+resWarn+","+resFail+"\n")
		return

# MAIN
if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='FASTQ Report')
	parser.add_argument('-n', '--name',help='Sample name',type=str, required=True)
	parser.add_argument('-R1', '--Read1',help='Fastq R1',type=str, required=True)
	parser.add_argument('-T1', '--Treads1',help='Fastq R1',type=str)
	args = parser.parse_args()

	if args.Treads1 != None:
		sample = Sample(args.name,args.Read1,args.Treads1)
		# RAW CHECK
		sample.raw_check()
		sample.raw_esito()
		
		# TRIMMED CHECK
		sample.trimm_check()
		sample.trimm_esito()
		
		# CREATE REPORT
		sample.makeReport()
	else:
		sample = Sample(args.name,args.Read1)
		# RAW CHECK
		sample.raw_check()
		sample.raw_esito()

		# CREATE REPORT
		sample.makeReport()
