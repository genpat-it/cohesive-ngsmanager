#!/usr/bin/env python3

import os
import sys


def process(ds, input_file, output_file):
    if os.path.isfile(input_file):
        with open(input_file, 'r') as input_file, open(output_file, 'w') as result_file:
            for line in input_file.readlines():
                split = line.strip().split("\t")
                cleaned = split[:-2]
                if ds in cleaned[0]:
                    for i in range(1, len(cleaned)):
                        if "+" in cleaned[i]:
                            cleaned[i] = cleaned[i].replace("+", "")
                        elif "-" in cleaned[i] or cleaned[i] == "N" or cleaned[i] == "0" or cleaned[i] == "0?":
                            cleaned[i] = "LNF"
                    result_file.write("\t".join(str(x) for x in cleaned) + "\n")
                else:
                    result_file.write("\t".join(str(x) for x in cleaned) + "\n")
    else:
        sys.stderr.write('could not open file: ' + input_file)
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.stderr.write('Usage: python cleanup.py ds inputFile outputFile \n')
        sys.exit(1)
    else:
        process(sys.argv[1], sys.argv[2], sys.argv[3])
