#!/usr/bin/env python3

import argparse


def select_coverage_file(coverage_files, threshold, output):
    selected_coverage = []
    for coverage_file in coverage_files:
        lines = coverage_file.readlines()
        if len(lines) != 4:
            exit(1)
        if float(lines[1].strip()) >= threshold:
            selected_coverage.append(lines)
    with open(output, 'w') as output_file:
        if len(selected_coverage) == 0:
            with open('errors.log', 'w') as errors_log:
                errors_log.write('no lineage could be found')
                output_file.write("lineage,reference,hcov\n")
                output_file.write("NA,-,-")    
        elif len(selected_coverage) > 1:
            with open('errors.log', 'w') as errors_log:
                errors_log.write('The sample cannot be associated to more than one lineage')
                output_file.write("lineage,reference,hcov\n")
                output_file.write("ND,-,-")    
        else:
                lines = selected_coverage[0]            
                output_file.write("lineage,reference,hcov\n")
                output_file.write("%s,%s,%s" % (lines[3].strip(), lines[2].strip(), lines[1].strip()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_list', dest='coverage_files', required=True,
                        help='coverage_files', nargs='+',
                        type=argparse.FileType('r'))
    parser.add_argument('--threshold', dest='threshold', required=True,
                        help='minimum hcov',
                        type=float)
    parser.add_argument('--output', dest='output', required=True,
                        help='output file name',
                        type=str)
    args = parser.parse_args()

    select_coverage_file(args.coverage_files, args.threshold, args.output)
