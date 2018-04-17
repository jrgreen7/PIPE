#!/usr/bin/python
#convert_names_to_IDs.py input_file protein_index_file output_file

import sys

if len(sys.argv) != 4:
    print 'USAGE: convert_names_to_IDs.py input_file protein_index_file output_file'
    print '\nconverts input_file from pairs of protein IDs to pairs of protein'
    print 'names (as defined by protein_index_file), written to output_file\n'
    sys.exit(0)

index_file = open(sys.argv[2], 'r')

proteins = dict([])

index_file.readline()

i = 0;
for line in index_file:
	proteins[line.split()[0].strip()] = i;
	i = i+1
	#proteins.append(line.split()[0].strip())

index_file.close()

input_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[3], 'w')

for line in input_file:
	output_file.write(str(proteins[line.split()[0].strip()]))
	output_file.write('\t')
	output_file.write(str(proteins[line.split()[1].strip()]))
	output_file.write('\n')

input_file.close()
output_file.close()
