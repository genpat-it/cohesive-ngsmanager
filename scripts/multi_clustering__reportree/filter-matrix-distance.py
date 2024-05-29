#!/usr/bin/env python3

import os
import sys


def distance_report(matrix_file, thresholds, samples, output_file):
    if os.path.isfile(matrix_file):
        with open(output_file, 'w') as res:
            for threshold in thresholds.split(','):            
                res.write("Distance report, threshold: {}\n\n".format(threshold))
                for sampleraw in samples.split(','):
                    sample = sampleraw.strip()
                    header = ''
                    with open(matrix_file, 'r') as dists:
                        for line in dists:
                            content_written = 0
                            row = line.strip().split('\t')
                            if header == '':
                                header = line.strip().split('\t')
                            elif row[0] == sample:
                                for col in range(len(header)):
                                    if col == 0:
                                        continue
                                    if header[col] == sample:
                                        continue
                                    if int(row[col]) <= int(threshold):
                                        if content_written == 0:
                                            res.write("Sample: {}\n".format(sample))
                                            content_written = 1
                                        res.write("{}\t{}\n".format(header[col], row[col]))
                                if content_written == 1:
                                    res.write("---\n")
                res.write("====\n\n")
    else:
        sys.stderr.write('could not open input files')
        sys.exit(1)
    with open(output_file, 'r') as res:
        print(res.read())


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.stderr.write(
            "Usage: python %s <cgmlst_dist_result> <threshold> <samples_of_interest_file> <output_file>" % (
                sys.argv[0]))
        sys.exit(1)
    else:
        distance_report(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
