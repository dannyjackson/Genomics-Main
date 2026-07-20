# Converting pyrho per-base per-generation recombination rate to CM/Mb and adding CM column

# This script first adds chrom_id column first to all rmaps
# Then it calculates the CM units using a file mapping some conversion rates for specific scaffolds. This assumes you have some external reference of recommended reference rates.
# These outputs can be easily concatenated together in a terminal (as in A2.5.1). These outputs are designed to be used in flexsweep with columns: [CHR_ID, Interval_start, Interval_end, CM/Mb, CM_distance]

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--cm_conversion_file", type=str, help="comma separated file mapping bp to cm conversion rates to specific scaffolds")
parser.add_argument("--pop", type=str, help="pop name prefix of rmaps. Assumes filenames are formatted as 'POP_SCAFFOLD.rmap'")
parser.add_argument("--indir", type=str, help="directory path of rmaps to parse")
args = parser.parse_args()

convert_file = pd.read_csv(args.cm_conversion_file)

outpath = os.path.join(args.indir, 'pyrho_cm_converted')
if not os.path.exists(outpath): os.mkdir(outpath)

# Generate converted files
for row in convert_file.index:
    scaffold = convert_file.loc[row, 'scaffold']
    convert_factor = convert_file.loc[row, 'cm_conversion_rate']
    
    print('Converting: ' + scaffold)
    rmap_path = os.path.join(args.indir, f'{args.pop}_{scaffold}.rmap')

    data = []
    with open(rmap_path, 'r') as file:
        for line in file:
            start, end, rate = line.strip().split()
            data.append([scaffold, int(start), int(end), float(rate)]) # adding chrom_id here
        
    rmap_df = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'bp_rate'])
    rmap_df['cm/mb'] = rmap_df['bp_rate'] * convert_factor
    rmap_df['cm'] = (rmap_df['start'] + rmap_df['end']) * rmap_df['cm/mb']
    
    rmap_df.drop(columns=['bp_rate'], inplace=True)
    rmap_df.to_csv(os.path.join(outpath, f'{args.pop}_{scaffold}_cm.rmap'), sep='\t', header=None, index=0)
    print(scaffold + ' converted pyrho rmap saved to: ' + os.path.join(outpath, f'{args.pop}_{scaffold}_cm.rmap'))