# Helper script to convert popsize CSV outputs from SMC++ to demes YAML format.
# Not a one-size-fits-all script. Edit as needed. This is just the basic idea

import os
import demes
#import demesdraw
import pandas as pd

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-c", "--csv_file", type=str, help='Path to CSV file output from SMC++ plotting')
parser.add_argument("-t", "--time_units", default='years', choices=['years', 'generations'], help="Specify time units. Either 'years', 'generations'. (Default to years)")
parser.add_argument("-d", "--description", type=str, default="popsizes from smc++", help='Add a description to the demes file')
parser.add_argument("-g", "--gen_time", type=int, help='Specify generation time. Be sure to set to 1 if using generations as time units')
parser.add_argument("-o", "--outfilename", type=str, help='Specify output yaml file name. File will output to same directory as scmc csv file')
args = parser.parse_args()

ne_df = pd.read_csv(args.csv_file)

b = demes.Builder(
    description=args.description,
    time_units=args.time_units,
    generation_time=args.gen_time,
)

for pop in set(ne_df['label']):
    epoch_lst = []
    for row in ne_df.index:
        if ne_df.loc[row,'label'] == pop:
            epoch_lst.insert(0, dict(start_size=ne_df.loc[row, 'y'], end_time=ne_df.loc[row, 'x'])) # Build list this way since demes wants ancient times first (reverse of SMC++)
    b.add_deme(pop, epochs=epoch_lst)

demes_graph = b.resolve()

#demesdraw.tubes(demes_graph)

outdir = os.path.split(args.csv_file)[0]
demes.dump(demes_graph, os.path.join(outdir, f'{args.outfilename}.yaml'))