# Reformatting each line of the .ped file as such (mainly need to just clean up the sample/chromosome ids)
# lamich PL15 0 0 0 -9 0 0 G G C T G G 0 0...     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      1 lamichPL15 0 0 0 -9 0 0 G G C T G G 0 0...

import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-p", "--ped", type=str, help='The path to the ped file being cleaned')
args = parser.parse_args()

with open(args.ped, 'r') as f:
    file = f.readlines()

cleaned_str_lst = []

for line in file:
    indv = line.split()
    if indv[0] != indv[1]:
        name = indv[0] + indv[1]
    else:
        name = indv[0]
    indv[0] = '1'
    indv[1] = name
    cleaned_str_lst.append(' '.join(indv)+'\n')
    #print(' '.join(indv[:50]))

with open(args.ped, 'w') as f:
    f.writelines(cleaned_str_lst)