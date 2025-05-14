#!/usr/bin/env python3

from Bio import SeqIO
from numpy import average
import os
import sys
import argparse
from multiprocessing import Pool
import functools

def decode(c):
    return ord(c) - 33	

def length_and_quality(title, sequence, quality):
    decoded_quality = list(map(decode, quality))
    qual = average(decoded_quality)
    q30_bases = sum(q >= 30 for q in decoded_quality)
    return [len(sequence), qual, q30_bases]

def safe_open(file, mode='rt'):
    # safe open file
    if file.endswith('.gz'):
        import gzip
        return gzip.open(file, mode)
    else:
        return open(file, mode)


def readsInfo(fastq):
    res = {}
    lengthReads = []
    total_bases = 0
    q30_bases = 0
    avgQual = []
    try:
        pool_size = min(os.cpu_count(), 8)

        with safe_open(fastq) as handle, Pool(pool_size) as pool:
            result = pool.starmap(length_and_quality, SeqIO.QualityIO.FastqGeneralIterator(handle))
        
        for (seqlen, qual, q30) in result:
            lengthReads.append(seqlen)
            avgQual.append(qual)
            total_bases += seqlen
            q30_bases += q30
            
        
        # lengthReads = list(map(lambda a: a[0], result)) 
        # avgQual = list(map(lambda a: a[1], result)) 
        # readsQual30 = len(list(filter(lambda x: x > 30, avgQual)))

    except Exception as e: 
        print(e)
        print("Errore nel file " + fastq)
        exit(1)

    if len(lengthReads) != 0:
        q30Perc = round(q30_bases / float(total_bases) * 100, 2)
        Mbases = round(float(sum(lengthReads))/1000000,2)
        res.update({"countReads":len(lengthReads)})
        res.update({"minLengthReads":min(lengthReads)})
        res.update({"maxLengthReads":max(lengthReads)})
        res.update({"avgLengthReads":average(lengthReads)})
        res.update({"minQual":min(avgQual)})
        res.update({"maxQual":max(avgQual)})
        res.update({"avgQual":average(avgQual)})
        res.update({"readsQ30": q30Perc})
        res.update({"Mbases":Mbases})
    else:
        res.update({"countReads":0})
        res.update({"minLengthReads":0})
        res.update({"maxLengthReads":0})
        res.update({"avgLengthReads":0})
        res.update({"minQual":0})
        res.update({"maxQual":0})
        res.update({"avgQual":0})
        res.update({"readsQ30":0})
        res.update({"Mbases":0})

    return res


# Class Definition
class SampleRaw:

    def __init__(self, read1, read2):
        # pool = mp.Pool(processes=2)
        print("Raw reads:")
        # fastqInfo = pool.map(readsInfo,[read1,read2])
        # self.R1 = fastqInfo[0]
        # self.R2 = fastqInfo[1]
        print("Processing R1...")
        self.R1 = readsInfo(read1)
        print("Processing R2...")
        self.R2 = readsInfo(read2)
        self.warnings = []
        self.fails = []

    def get_R1(self):
        return self.R1

    def get_R2(self):
        return self.R2

    def get_warnings(self):
        return self.warnings

    def get_fails(self):
        return self.fails

    def getNR(self):
        totReads = self.R1.get("countReads") + self.R2.get("countReads")
        return totReads

    def get_q30Reads(self):
        totReads = self.R1.get("countReads") + self.R2.get("countReads")
        totQ30 = self.R1.get("readsQ30") + self.R2.get("readsQ30")
        if float(totReads) == 0:
            return 0
        return round(float(totQ30)/2,2)

    def get_stats(self):
        minLength = min(self.R1.get("minLengthReads"), self.R2.get("minLengthReads"))
        maxLength = min(self.R1.get("maxLengthReads"), self.R2.get("maxLengthReads"))
        avgLength = round((self.R1.get("avgLengthReads") + self.R2.get("avgLengthReads")) / 2, 2)
        minQual = min(self.R1.get("minQual"), self.R2.get("minQual"))
        maxQual = min(self.R1.get("maxQual"), self.R2.get("maxQual"))
        avgQual = round((self.R1.get("avgQual") + self.R2.get("avgQual")) / 2, 2)
        totMbases = self.R1.get("Mbases") + self.R2.get("Mbases")
        return [self.getNR(), totMbases, minLength, maxLength, avgLength, self.get_q30Reads(), minQual, maxQual,
                avgQual]

    def checkPair(self):
        if self.R1.get("countReads") != self.R2.get("countReads"):
            self.fails.append("FAIL_raw_pairNR")

    def checkQ30(self):
        q30R = round((self.R1.get("avgQual") + self.R2.get("avgQual")) / 2, 2)
        if float(q30R) >= 25:
            return
        elif float(q30R) >= 5 and float(q30R) < 25:
            self.warnings.append("WARN_rawQ30")
        else:
            self.fails.append("FAIL_rawQ30")

    def checkQavg(self):
        avgQualR = round((self.R1.get("avgQual") + self.R2.get("avgQual")) / 2, 2)
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

    def __init__(self, read1, read2, unp):
        # pool = mp.Pool(processes=3)
        print("Trimmed reads:")
        # fastqInfo = pool.map(readsInfo,[read1,read2,unp])
        # self.R1 = fastqInfo[0]
        # self.R2 = fastqInfo[1]
        # self.U = fastqInfo[2]
        print("Processing R1...")
        self.R1 = readsInfo(read1)
        print("Processing R2...")
        self.R2 = readsInfo(read2)
        print("Processing Unpaired...")
        self.U = readsInfo(unp)
        self.warnings = []
        self.fails = []

    def get_R1(self):
        return self.R1

    def get_R2(self):
        return self.R2

    def get_unpaired(self):
        return self.U

    def get_warnings(self):
        return self.warnings

    def get_fails(self):
        return self.fails

    def getNR_total(self):
        totReads = self.R1.get("countReads") + self.R2.get("countReads") + self.U.get("countReads")
        return totReads

    def getNR(self):
        totReads = self.R1.get("countReads") + self.R2.get("countReads")
        return totReads

    def getNR_unpaired(self):
        return self.U.get("countReads")

    def get_q30Reads(self):
        totReads = self.R1.get("countReads") + self.R2.get("countReads")
        totQ30 = self.R1.get("readsQ30") + self.R2.get("readsQ30")
        if float(totReads) == 0:
            return 0
        return round(float(totQ30)/2,2)

    def get_stats(self):
        minLength = min(self.R1.get("minLengthReads"), self.R2.get("minLengthReads"))
        maxLength = min(self.R1.get("maxLengthReads"), self.R2.get("maxLengthReads"))
        avgLength = round((self.R1.get("avgLengthReads") + self.R2.get("avgLengthReads")) / 2, 2)
        minQual = min(self.R1.get("minQual"), self.R2.get("minQual"))
        maxQual = min(self.R1.get("maxQual"), self.R2.get("maxQual"))
        avgQual = round((self.R1.get("avgQual") + self.R2.get("avgQual")) / 2, 2)
        totMbases = self.R1.get("Mbases") + self.R2.get("Mbases")
        return [self.getNR_total(), self.getNR_unpaired(), self.getNR(), totMbases, minLength, maxLength, avgLength,
                self.get_q30Reads(), minQual, maxQual, avgQual]

    def checkPair(self):
        if self.R1.get("countReads") != self.R2.get("countReads"):
            self.fails.append("FAIL_trim_pairNR")

    def checkQ30(self):
        q30R = round((self.R1.get("avgQual") + self.R2.get("avgQual")) / 2, 2)
        if float(q30R) >= 35:
            return
        elif float(q30R) >= 10 and float(q30R) < 35:
            self.warnings.append("WARN_trimQ30")
        else:
            self.fails.append("FAIL_trimQ30")

    def checkQavg(self):
        avgQualR = round((self.R1.get("avgQual") + self.R2.get("avgQual")) / 2, 2)
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

    def __init__(self, sampleName, R1, R2, T1, T2, Unp):
        self.name = sampleName
        self.Raw = SampleRaw(R1, R2)
        self.Trimmed = SampleTrimmed(T1, T2, Unp)
        self.warnings = []
        self.fails = []

    def get_name(self):
        return self.name

    def raw_check(self):
        self.Raw.checkPair()
        self.Raw.checkQ30()
        self.Raw.checkQavg()
        self.Raw.checkNR_Bact()
        self.Raw.checkNR_Vir()
        with open(self.get_name() + "_SRC_raw.csv", 'w') as checkOut:
            results = ",".join(str(x) for x in self.Raw.get_stats())
            checkOut.write(
                "#Sample_name,Total_reads,Total_Mbases,Min_length,Max_length,Mean_length,Q30_reads,Min_avgQual,Max_avgQual,Mean_avgQual\n")
            checkOut.write(self.get_name() + "," + results + "\n")
        return

    def raw_esito(self):
        resWarn = ":".join(self.Raw.get_warnings())
        resFail = ":".join(self.Raw.get_fails())
        with open(self.get_name() + "_SRC_raw.check", 'w') as resOut:
            resOut.write("#Sample_name,WARNING,FAIL\n")
            resOut.write(self.get_name() + "," + resWarn + "," + resFail + "\n")
        return

    def trimm_check(self):
        self.Trimmed.checkPair()
        self.Trimmed.checkQ30()
        self.Trimmed.checkQavg()
        self.Trimmed.checkNR_Bact()
        self.Trimmed.checkNR_Vir()
        with open(self.get_name() + "_SRC_treads.csv", 'w') as checkOut:
            results = ",".join(str(x) for x in self.Trimmed.get_stats())
            checkOut.write(
                "#Sample_name,Total_reads,Unpaired_reads,Paired_reads,Paired_Mbases,Min_PairedLength,Max_PairedLength,Mean_PairedLength,Q30_PairedReads,Min_PairedAvgQual,Max_PairedAvgQual,Mean_PairedAvgQual\n")
            checkOut.write(self.get_name() + "," + results + "\n")
        return

    def trimm_esito(self):
        resWarn = ":".join(self.Trimmed.get_warnings())
        resFail = ":".join(self.Trimmed.get_fails())
        with open(self.get_name() + "_SRC_treads.check", 'w') as resOut:
            resOut.write("#Sample_name,WARNING,FAIL\n")
            resOut.write(self.get_name() + "," + resWarn + "," + resFail + "\n")
        return

    def trimDiscarded(self):
        discarded = self.Raw.getNR() - self.Trimmed.getNR_total()
        return discarded

    def makeReport(self):
        with open(self.get_name() + "_readsCheck.csv", 'w') as checkOut:
            resRaw = self.Raw.get_stats()
            resTrim = self.Trimmed.get_stats()
            resWarn = ":".join(self.Raw.get_warnings() + self.Trimmed.get_warnings())
            resFail = ":".join(self.Raw.get_fails() + self.Trimmed.get_fails())
            results = ",".join(str(x) for x in
                               [resRaw[0], resRaw[1], resRaw[4], resRaw[5], resRaw[8], resTrim[0], self.trimDiscarded(),
                                resTrim[1], resTrim[2], resTrim[3], resTrim[6], resTrim[7], resTrim[10], 0])
            checkOut.write(
                "#Sample_name,Total_rawReads,Total_rawMbases,Mean_rawLength,Q30_rawReads,Mean_rawAvgQual,Total_trimReads,trimDiscarded,trimUnpaired,trimPaired,trimPairedMbases,Mean_trimPairedLength,Q30_trimPairedReads,Mean_trimPairedAvgQual,Overlap_trimPaired,WARNING,FAIL\n")
            checkOut.write(self.get_name() + "," + results + "," + resWarn + "," + resFail + "\n")
        return


# MAIN
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FASTQ Report')
    parser.add_argument('-n', '--name', help='Sample name', type=str, required=True)
    parser.add_argument('-R1', '--Read1', help='Fastq R1', type=str, required=True)
    parser.add_argument('-R2', '--Read2', help='Fastq R2', type=str, required=True)
    parser.add_argument('-T1', '--Treads1', help='Fastq trimmed R1', type=str, required=True)
    parser.add_argument('-T2', '--Treads2', help='Fastq trimmed R2', type=str, required=True)
    parser.add_argument('-U', '--Unpaired', help='Fastq unpaired', type=str, required=True)
    args = parser.parse_args()

    sample = Sample(args.name, args.Read1, args.Read2, args.Treads1, args.Treads2, args.Unpaired)

    # RAW CHECK
    print("Raw reads CHECK")
    sample.raw_check()
    print("Raw reads REPORT")
    sample.raw_esito()

    # TRIMMED CHECK
    print("Trimmed reads CHECK")
    sample.trimm_check()
    print("Trimmed reads REPORT")
    sample.trimm_esito()

    # CREATE REPORT
    print("Create final report")
    sample.makeReport()
    print("DONE")