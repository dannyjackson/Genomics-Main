# Converting pyrho per-base per-generation recombination rate to CM/Mb and adding CM column

# This script first adds chrom_id column first to all rmaps
# Then it calculates the CM units
# These outputs can be easily concatenated together in a terminal (as in A2.5.1). These outputs are designed to be used in flexsweep with columns: [CHR_ID, Interval_start, Interval_end, CM/Mb, CM_distance]

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--cm_conversion", type=float, help="CM conversion rate", default=1e8)
parser.add_argument("--scaffolds", type=str, help="scaffolds list file")
parser.add_argument("--pop", type=str, help="pop name prefix of rmaps. Assumes filenames are formatted as 'POP_SCAFFOLD.rmap'")
parser.add_argument("--indir", type=str, help="directory path of rmaps to parse")
args = parser.parse_args()

with open(args.scaffolds, 'r') as f:
    scaffolds = [line.strip() for line in f.readlines()]

outpath = os.path.join(args.indir, 'pyrho_cm_converted')
if not os.path.exists(outpath): os.mkdir(outpath)

# Generate converted files
for s in scaffolds:
    print('Converting: ' + s)
    rmap_path = os.path.join(args.indir, f'{args.pop}_{s}.rmap')

    data = []
    with open(rmap_path, 'r') as file:
        for line in file:
            start, end, rate = line.strip().split()
            data.append([s, int(start), int(end), float(rate)]) # adding chrom_id here
        
    rmap_df = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'bp_rate'])
    rmap_df['cm/mb'] = rmap_df['bp_rate'] * args.cm_conversion
    rmap_df['cm'] = (rmap_df['start'] + rmap_df['end']) * rmap_df['cm/mb']
    
    rmap_df.drop(columns=['bp_rate'], inplace=True)
    rmap_df.to_csv(os.path.join(outpath, f'{args.pop}_merged.rmap'), sep='\t', header=None, index=0)
    print('Merged pyrho rmap saved to: ' + os.path.join(outpath, f'{args.pop}_merged.rmap'))