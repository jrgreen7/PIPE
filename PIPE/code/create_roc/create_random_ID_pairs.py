#!/usr/bin/python
#create_random_pairs.py known_ID_pairs_file total_num_proteins num_pairs_to_produce output_file

import sys
import random

if len(sys.argv) != 5:
    print 'USAGE: create_random_pairs.py known_ID_pairs_file total_num_proteins num_pairs_to_produce output_file'
    print '\ncreates num_pairs_to_produce random protein ID pairs (IDs range to total_num_proteins). Script'
    print  'ensures pairs generated don\'t occur in known_ID_pairs_file. Output pairs written to output_file\n'
    sys.exit(0)


known_pairs = set([])

kp_file = open(sys.argv[1], 'r')

for line in kp_file:
    cols = line.strip().split('\t')
    if len(cols) == 2:
        known_pairs.add(cols[0] + '\t' + cols[1])

kp_file.close()

proteins = []

for i in range(int(sys.argv[2])):
    proteins.append(str(i))

r1 = random.SystemRandom()
r2 = random.SystemRandom()
num_pairs = int(sys.argv[3])

pairs = set([])

output_file = open(sys.argv[4], 'w')

while len(pairs) < num_pairs:
    p1 = r1.randint(0, len(proteins)-1)
    p2 = r2.randint(0, len(proteins)-1)
	
    pair = str(proteins[p1]) + '\t' + str(proteins[p2])
    pair2 = str(proteins[p2]) + '\t' + str(proteins[p1])
	
    if pair not in known_pairs and pair2 not in known_pairs:
        if pair not in pairs and pair2 not in pairs:
            if len(pairs) == 0:
                output_file.write(str(num_pairs) + '\n')
            pairs.add(pair)
            output_file.write(pair + '\n')
output_file.close()
		
