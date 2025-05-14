#!/usr/bin/env python3

import argparse
import pandas as pd

def find_closest_samples(matrix_file, samples, closest_num, output_file, include_input):
    # Load the distance matrix
    df = pd.read_csv(matrix_file, sep="\t", index_col=0)

    # Convert all values to numeric (strict conversion, will fail if any non-numeric data is present)
    df = df.apply(pd.to_numeric, errors='raise')

    # Trim whitespace and split sample codes
    target_samples = {sample.strip() for sample in samples.split(",")}

    # Set to store unique closest samples
    closest_samples = set()

    for sample in target_samples:
        if sample not in df.index:
            raise ValueError(f"Error: Sample '{sample}' not found in the matrix.")

        # Get distances, sort numerically, and select closest
        closest = df.loc[sample].sort_values().iloc[1:closest_num + 1].index
        closest_samples.update(closest)

    # If --include_input is set, add original input samples to the output
    if include_input:
        closest_samples.update(target_samples)

    # Save unique closest samples
    with open(output_file, "w") as f:
        f.write("\n".join(closest_samples))

    print(f"Saved {len(closest_samples)} unique closest samples to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find closest samples in a distance matrix.")
    parser.add_argument("--matrix", required=True, help="Path to the distance matrix file")
    parser.add_argument("--samples", required=True, help="Comma-separated list of sample codes")
    parser.add_argument("--closest_num", type=int, required=True, help="Number of closest samples per input sample")
    parser.add_argument("--output", required=True, help="Output file for closest samples")
    parser.add_argument("--include_input", action="store_true", help="Include input samples in the output file")

    args = parser.parse_args()
    find_closest_samples(args.matrix, args.samples, args.closest_num, args.output, args.include_input)
