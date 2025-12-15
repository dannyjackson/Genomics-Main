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
    clean_lst = []
    splt = line.split()
    if splt[0] not in parsed_scaffolds:
        scaffold_num += 1
        parsed_scaffolds.append(splt[0])
    snp_name = '.'.join([splt[0], splt[3]])
    genetic_loc_cm = splt[2]
    genetic_loc_snp = splt[3]
    clean_lst.append(str(scaffold_num))
    clean_lst.append(snp_name)
    clean_lst.append(genetic_loc_cm)
    clean_lst.append(genetic_loc_snp)
    #print('\t'.join(clean_lst))
    cleaned_map.append('\t'.join(clean_lst)+'\n')

with open(args.map, 'w') as f:
    mapfile = f.writelines(cleaned_map)