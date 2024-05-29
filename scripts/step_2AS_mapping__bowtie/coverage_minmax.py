#!/usr/bin/env python3
import sys


def create_coverage_report(sample, ds, samtools_view_output, samtools_depth_output, import_file_result):
    v_file = open(samtools_view_output, "r")
    d_file = open(samtools_depth_output, "r")
    num_mapped_reads = v_file.readlines()[0].strip()
    # get [min, max] vcov   
    if num_mapped_reads == '0':
        min_vcov = 0
        max_vcov = 0
    else:
        min_vcov = sys.maxsize
        max_vcov = -1
        for line in d_file.readlines():
            split_line = line.strip().split("\t")
            if len(split_line) > 2:
                value = int(split_line[2])
                if value > max_vcov:
                    max_vcov = value
                if value < min_vcov:
                    min_vcov = value    
        if min_vcov == sys.maxsize:
            min_vcov = 0
        if max_vcov == -1:
            min_vcov = 0
    result = open(import_file_result, 'w')
    result.write("CMP_ID,SAMPLE_DS,NUM_MAPPED_READS,MIN_VCOV,MAX_VCOV\n")
    result.write("{},{},{},{},{}\n".format(sample, ds.replace("DS", ""), num_mapped_reads, min_vcov, max_vcov))
    result.close()
    v_file.close()
    d_file.close()


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.stderr.write(
            "Usage: %s <sample> <ds> <samtools_view_out> <samtools_depth_out> <import_file_result>" % (sys.argv[0]))
        sys.exit(1)
    else:
        create_coverage_report(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
