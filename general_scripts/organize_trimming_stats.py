# A quick script to grab the trimming stats from the slurm output in A0.1_trimming.sh
# Better and less storage intensive than using -trimlog argument
# We're just searching for some specific strings outputted by trimmomatic that are always the same and using those to parse the read stats
# Places trimming stats csv file into working directory

slurm_outs = 'path/to/trimming/slurm/out/files'
sample_ids = 'path/to/sampleid/list/'

import os
import pandas as pd

trim_df = pd.DataFrame(columns=['Sample','both_survive', 'forward_survive', 'reverse_survive', 'dropped', 'inputted_pairs'])
with open(sample_ids) as f:
    ids = f.readlines()
ids = [id.strip() for id in ids]
trim_df['Sample'] = ids

for file in os.listdir(slurm_outs):
    with open(os.path.join(slurm_outs, file)) as f:
        outfile = f.readlines()
    for line in outfile:
        if 'Input Read Pairs:' in line:
            stats = line.strip()
    splt_stats = stats.split(')')
    dropped = float(splt_stats[3].split(':')[1].split('(')[0].strip())
    r_survive = float(splt_stats[2].split(':')[1].split('(')[0].strip())
    f_survive = float(splt_stats[1].split(':')[1].split('(')[0].strip())
    input_pairs = float(splt_stats[0].split('Both Surviving:')[0].split(':')[1].split('(')[0].strip())
    both_survive = float(splt_stats[0].split('Both Surviving:')[1].split('(')[0].strip())
    sample_num = file.split('_')[-1].strip('.out')
    info_lst = [sample_num, both_survive, f_survive, r_survive, dropped, input_pairs]

    idx = int(info_lst[0])-1
    dct = {'Sample':trim_df.loc[idx,'Sample'] ,'sample_num':info_lst[0], 'both_survive':info_lst[1], 'forward_survive':info_lst[2], 
           'reverse_survive':info_lst[3], 'dropped':info_lst[4], 'inputted_pairs':info_lst[5]}
    trim_df.loc[idx] = dct


trim_df.to_csv('trimming_stats.csv')
