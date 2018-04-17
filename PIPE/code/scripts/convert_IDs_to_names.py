#!/usr/bin/python
#convert_names_to_IDs.py input_file protein_index_file output_file

import sys

if len(sys.argv) != 4:
    print 'USAGE: convert_names_to_IDs.py input_file protein_index_file output_file'
    print '\nconverts input_file from pairs of protein names to pairs of'
    print 'protein IDs (as defined by protein_index_file), written to output_file\n'
    sys.exit(0)

index_file = open(sys.argv[2], 'r')

proteins = []

index_file.readline()

for line in index_file:
	proteins.append(line.split()[0].strip())

index_file.close()

input_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[3], 'w')

for line in input_file:
	cols = line.split('\t')
	output = ''

	for col in cols:
		output = output + proteins[int(col)] + '\t'
	output = output.strip()
	output_file.write(output + '\n')

input_file.close()
output_file.close()
