# Small helper script to reformat the chromosome conversion file built in base_setup.sh

import pandas as pd
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help='The path to the tab-separated chrom name map file generated from the NCBI dataset/dataformat commands')
parser.add_argument("-o", "--output", type=str, help='The path for the chrom_conversion file')
parser.add_argument("-e", "--exclusions", type=str, help='Strings to exclude, comma seperated, no spaces (usually portions of NCBI accession codes for incomplete scaffolds or sex chromosomes)')
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', skiprows=1, header=None)

# Split string of scaffold/chromosome codes we want into a list
exclusion_lst = args.exclusions.split(',')

for row in df.index:
    if any(code in df.loc[row,1] for code in exclusion_lst):
        df.drop(row, inplace=True)

df.set_index(0, inplace=True)
df.to_csv(args.output, header=None)