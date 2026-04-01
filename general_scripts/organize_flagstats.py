import json
import pandas as pd
import os

# This script is meant aid in organizing the flagstats json output (from A0.2_alignandsort.sh) into a single datafile. 

qc_type=''
fname_suffix = ''
fstat_name_dct = {'total': 'total_reads','primary': 'primary','secondary': 'secondary','supplementary': 'supplementary','duplicates': 'duplicates','primary duplicates': 'primary_duplicates',
 'mapped': 'mapped','mapped %': 'percent_mapped','primary mapped': 'primary_mapped','primary mapped %': 'percent_primary_mapped','paired in sequencing': 'paired_in_seq',
 'read1': 'read1','read2': 'read2','properly paired': 'proper_pair','properly paired %': 'precent_proper_pair','with itself and mate mapped': 'itself_and_mate','singletons': 'singleton','singletons %': 'percent_singleton',
 'with mate mapped to a different chr': 'mate_mapped_diff_chr','with mate mapped to a different chr (mapQ >= 5)': 'mate_mapped_diff_chr_mapQ>5'}

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-f", "--fstats", type=str, nargs='+', default=['percent_mapped', 'percent_primary_mapped', 'percent_proper_pair'], choices=list(fstat_name_dct.values()), help='Specify what flagstats to output (Default to percentages of mapped, primary mapped, and proper pair info')
parser.add_argument("-t", "--stattype", default='pass', choices=['pass', 'fail', 'both'], help="Specify what stats type to output. Either 'pass', 'fail', or 'both'. (Default to pass)")
parser.add_argument("-a", "--all", action='store_true', help='Specify to output ALL stats')
parser.add_argument("-s", "--samples", type=str, help='Specify path to file with sample names to collect flagstats from. Optional for if only want a subset of samples to be gathered')
parser.add_argument("-d", "--directory", type=str, help="Specify directory containing where flagstats are. Should be the 'datafiles/sortedbamfiles/' directory if following Gen-Main pipeline")
args = parser.parse_args()

def filter_cols(df, opts):
    for col in df:
        if col not in opts:
            df.drop(col, axis=1)

def rename_cols(df):
    for col in df.columns:
        if 'passed' in col:
            df.rename(columns={col:f'passed.{col.split('.')[1]}'}, inplace=True)
        elif 'failed' in col: 
            df.rename(columns={col:f'failed.{col.split('.')[1]}'}, inplace=True)

# Decide stattype
if not args.all:
    match args.stattype:
        case 'pass':
            qc_type = "QC-passed reads"
            fname_suffix = '_passed'
            print('Grabbing passed read info')
        case 'fail':
            qc_type = "QC-failed reads"
            fname_suffix = '_failed'
            print('Grabbing failed read info')
        case 'both':
            qc_type = "both"
            print('Grabbing passed and failed read info')


# Decide sample list
if args.samples:
    print('User inputted sample file')
    sample_lst = args.samples
else:
    print('Assuming sample names by files in given directory')
    sample_lst = os.listdir(args.directory)

# Print flagstat options
if args.all:
    fname_suffix = '_all'
    print('Outputting ALL flagstats to datafile') 
else:
    print(f'Requested Data:{args.fstats}')


# Grab stats for each sample
fstat_df_lst = []
for sample in sample_lst:
    with open(os.path.join(args.directory, sample, f'{sample}_flagstat.json')) as j:
        sample_file = json.load(j)

    if args.all == True or qc_type == 'both':
        sample_stats = pd.json_normalize(sample_file)
        rename_cols(sample_stats)
        if qc_type == 'both': filter_cols(sample_stats, args.fstats)
    else:
        sample_stats = pd.json_normalize(sample_file[qc_type])
        sample_stats.rename(columns=fstat_name_dct, inplace=True)
        filter_cols(sample_stats, args.fstats)
    
    fstat_df_lst.append(sample_stats)

fstat_df = pd.concat(fstat_df_lst, ignore_index=True)
fstat_df.to_csv(os.path.join(args.directory, 'flagstats', f'flagstats{fname_suffix}.csv'))

print('Finished grabbing flagstats')