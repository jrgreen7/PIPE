#!/usr/bin/python
#check_sequences_and_pairs.py sequence_file pairs_file

import sys

if len(sys.argv) != 3:
    print 'USAGE: check_sequences_and_pairs.py sequence_file pairs_file'
    print '\nensures that the sequence file and pairs file formatted correctly,'
    print 'each sequence contains valid amino acids and each pair contains'
    print 'proteins found in the sequence file\n'
    sys.exit(0)

AAs = ['A','R', 'N','D', 'C', 'Q','E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S','T', 'W', 'Y', 'V', 'B', 'Z', 'X', 'U']
longest_protein = 0;
proteins = set([])

sequence_file = open(sys.argv[1], 'r')
    
for line in sequence_file:
    cols = line.strip().split('\t')

    if len(cols) != 2:
        print 'ERROR: Sequence file not formatted properly'
        sys.exit(0)

    proteins.add(cols[0])

    for char in cols[1]:
        if char not in AAs and char != '\n':
            print 'ERROR: Bad sequence'
            sys.exit(0)

    if len(cols[1]) > longest_protein:
        longest_protein = len(cols[1])

print 'max protein length: ' + str(longest_protein)
sequence_file.close()

num_known_pairs = 0
pair_file = open(sys.argv[2], 'r')
for line in pair_file:
    cols = line.strip().split('\t')

    if len(cols) != 2:
        print 'ERROR: Protein pairs file not formatted properly'
        sys.exit(0)

    if cols[0] not in proteins:
        print 'ERROR: Protein ' + cols[0] + ' appears in protein pairs but not in sequence file'
        sys.exit(0)

    if cols[1] not in proteins:
        print 'ERROR: Protein ' + cols[1] + ' appears in protein pairs but not in sequence file'
        sys.exit(0)

    num_known_pairs = num_known_pairs + 1

print 'num known pairs: ' + str(num_known_pairs)
pair_file.close()
