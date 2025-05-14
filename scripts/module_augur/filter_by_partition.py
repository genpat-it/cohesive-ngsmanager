#!/usr/bin/env python3

import pandas as pd
import argparse


def extract_codes(tsv_file, partition, output_file):
    # Load TSV file into a DataFrame
    df = pd.read_csv(tsv_file, sep='\t')

    # Filter rows where partition matches the parameter
    filtered_df = df[df['partition'] == partition]

    # Extract and split codes from 'samples' column
    unique_codes = set()
    for samples in filtered_df['samples'].dropna():  # Drop NaN values
        unique_codes.update(samples.split(','))

    # Write unique codes to output file
    with open(output_file, 'w') as f:
        for code in sorted(unique_codes):
            f.write(code + '\n')

    # Print log message
    print(f"Selected {len(unique_codes)} samples from clusters related to SOI for partition: {partition}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract unique sample codes from a TSV file based on partition.")
    parser.add_argument("--partition_summary", required=True, help="Path to the TSV file.")
    parser.add_argument("--partition", required=True, help="Partition to filter.")
    parser.add_argument("--output_samples", required=True, help="Output file for extracted sample codes.")

    args = parser.parse_args()

    extract_codes(args.partition_summary, args.partition, args.output_samples)