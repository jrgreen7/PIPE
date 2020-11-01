#!/usr/bin/python
""" setup_agcan.py
Author: Kevin Dick
Date: 2018-12-11
---
Setup each of the seven AGCAN jobs.:
Recompiles their own specific version of PIPE
given their unique parameters.
NOTE: Can only be run when all database files generated.........
NOTE2: Proceeding with only Intras and doing Inters later...
"""
import os, sys
import glob

jobs = ['cvd']

base_dir = '/home/williamm/scratch/Deep-PIPE-Sites/PIPE4/'

for job_id in jobs:
    print 'Running Setup for Job: ' + job_id

    data_dir = base_dir + '/'
    code_dir = base_dir + '/code/'
    

    print 'Checking protein sequences...'
    AAs = ['A','R', 'N','D', 'C', 'Q','E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S','T', 'W', 'Y', 'V', 'B', 'Z', 'X', 'U']
    proteins = set([])
    longest_protein = 0;
    protein_index = {}

    sequence_file = open(data_dir + 'data/protein_sequences.txt', 'r')
    index = 0
    for line in sequence_file:
        cols = line.strip().split('\t')

        if len(cols) != 2:
            print line
            print 'ERROR: Sequence file not formatted properly'
            sys.exit(0)

        proteins.add(cols[0])
        protein_index[cols[0]] = index
        index += 1

        for char in cols[1]:
            if char not in AAs and char != '\n':
                print 'ERROR: Bad sequence'
		print 'Char: ' + char + ' | Index: ' + str(index)
                sys.exit(0)

        if len(cols[1]) > longest_protein:
            longest_protein = len(cols[1])

    sequence_file.close()


    print 'Checking protein pairs...'
    num_known_pairs = 0

    pair_file = open(data_dir + 'data/protein_pairs.txt', 'r')
    for line in pair_file:
        cols = line.strip().split('\t')

        if len(cols) != 2:
            print 'ERROR: Protein pairs file not formatted properly'
            print line
	    sys.exit(0)

        if cols[0] not in proteins:
            print 'ERROR: Protein ' + cols[0] + 'appears in protein pairs but not in sequence file'
            sys.exit(0)

        if cols[1] not in proteins:
            print 'ERROR: Protein ' + cols[1] + 'appears in protein pairs but not in sequence file'
            sys.exit(0)

        num_known_pairs = num_known_pairs + 1
    pair_file.close()

    ###################################
    ## Convert pairs file to indexes ##
    ###################################
    print 'Creating protein indexES file...'
    with open(data_dir + 'data/protein_pairs.txt', 'r') as in_pairs, open(data_dir + 'data/protein_pairs_indexes.txt','w') as out_pairs:
        out_pairs.write("%d\n" % num_known_pairs)
        for line in in_pairs:
            cols = line.strip().split('\t')
            out_pairs.write('%d\t%d\n' % (protein_index[cols[0]], protein_index[cols[1]]))

    ################################################ # THIS WAS ACCOMPLISHED PREVIOUSLY...
    ## Create the interaction graph file from the ##
    ## protein pairs and sequence files           ##
    ################################################
    print 'Creating protein index file...'
    os.system('cp ' + code_dir + 'genTab/convertPairs.pl ' + data_dir + 'data/')
    os.system(data_dir + 'data/convertPairs.pl')
    os.system('rm ' + data_dir + 'data/convertPairs.pl')

    ###################################################
    ## Determine all other necessary PIPE parameters ##
    ###################################################
    num_set_bits = 0
    while num_set_bits < len(proteins):
        num_set_bits = num_set_bits + 64

    max_neighbours = -10
    f = open(data_dir + 'data/protein_pairs_index.txt', 'r')
    for line in f:
        if len(line.strip().split()) - 2 > max_neighbours:
            max_neighbours = len(line.strip().split()) - 2
    f.close()

    max_db_file = 0
    db_files = glob.glob(data_dir + 'database/*')
    for db_file in db_files:
        if os.path.getsize(db_file) > max_db_file:
            max_db_file = os.path.getsize(db_file)

    #############################
    ## Write paramters to file ##
    #############################
    print 'Writing Job Parameters to File...'
    f = open(data_dir + 'parameters.txt', 'w')

    f.write('num proteins\t' + str(len(proteins)) + '\n')
    f.write('num known pairs\t' + str(num_known_pairs) + '\n')
    f.write('num pairs\t' + str((len(proteins) * (1+len(proteins)))/2) + '\n')
    f.write('max_db_file\t' + str(max_db_file) + '\n')
    f.write('num_set_bits\t' + str(num_set_bits) + '\n')
    f.write('max_protein_len\t' + str(longest_protein) + '\n')
    f.write('max_neighbours\t' + str(max_neighbours) + '\n')
    f.close()

    #############################################################################
    ## Update the mp-pipe2 code to reflect the new parameters and recompile it ##
    #############################################################################
    cmd = 'cp ' + code_dir + 'MP-PIPE2/mp-pipe2.c ' + code_dir + 'MP-PIPE2/mp-pipe2.c_bak'
    print cmd
    os.system(cmd)
    in_file  = open(code_dir + 'MP-PIPE2/mp-pipe2.c_bak', 'r')
    out_file = open(code_dir + 'MP-PIPE2/mp-pipe2.c', 'w')

    for line in in_file:
        if line.startswith('#define MAX_DB_FILE'):
            out_file.write('#define MAX_DB_FILE ' + str(max_db_file) + '\n')
        elif line.startswith('#define NUM_SET_BITS'):
            out_file.write('#define NUM_SET_BITS ' + str(num_set_bits) + '\n')
        elif line.startswith('#define MAX_PROTEIN_LEN'):
            out_file.write('#define MAX_PROTEIN_LEN ' + str(longest_protein) + '\n')
        elif line.startswith('#define MAX_NEIGHBOURS'):
            out_file.write('#define MAX_NEIGHBOURS ' + str(max_neighbours) + '\n')
        else:
            out_file.write(line)
    in_file.close()
    out_file.close()

    # RECOMPILE
    print 'Recompiling...'
    cmd = 'mpicc -O3 -fopenmp -Wall ' + code_dir + 'MP-PIPE2/mp-pipe2.c -m64 -lm -o ' + code_dir + 'MP-PIPE2/mp-pipe2'
    print cmd
    os.system(cmd)

print 'Execution Complete!~'
