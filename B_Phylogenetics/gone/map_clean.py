# Reformatting each line of the .map file as such (basically reordering columns)
#NC_044571.1	.	0	450     >>>>>>>>>>>>>>>>   1	NC_044571.1.450	0	450

import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-m", "--map", type=str, help='The path to the map file being cleaned')
args = parser.parse_args()

with open(args.map, 'r') as f:
    mapfile = f.readlines()

cleaned_map = []
parsed_scaffolds = []
scaffold_num = 0

for line in mapfile:
    snp_info = line.split()

    if snp_info[0] not in parsed_scaffolds:
        parsed_scaffolds.append(snp_info[0])
        scaffold_num += 1
    snp_name = '.'.join([snp_info[0], snp_info[3]])
    genetic_loc_cm = snp_info[2]
    genetic_loc_snp = snp_info[3]

    cleaned_map.append('\t'.join([str(scaffold_num), snp_name, genetic_loc_cm, genetic_loc_snp])+'\n')

with open(args.map, 'w') as f:
    mapfile = f.writelines(cleaned_map)

print('Completed map file cleanup')