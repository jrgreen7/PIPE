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

jobs = ['yeast']

base_dir = '/home/williamm/scratch/final-wrap-up/PIPE/PIPE4/' # do not include trailing slash

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
    ## backup convertPairs.pl
    os.system('cp ' + code_dir + 'genTab/convertPairs.pl ' + code_dir + 'genTab/convertPairs.pl_bak')
    in_file  = open(code_dir + 'genTab/convertPairs.pl_bak', 'r')
    out_file = open(code_dir + 'genTab/convertPairs.pl', 'w')
    for line in in_file:
        if line.startswith('$SEQ_FILE'):
            out_file.write('$SEQ_FILE=\"' + data_dir + 'data/protein_sequences.txt' + '\";\n')
        elif line.startswith('$STR_PAIR_FILE'):
            out_file.write('$STR_PAIR_FILE=\"' + data_dir + 'data/protein_pairs.txt' + '\";\n')
        elif line.startswith('$IDX_PAIR_FILE'):
            out_file.write('$IDX_PAIR_FILE=\"' + data_dir + 'data/protein_pairs_index.txt' + '\";\n')
        else:
            out_file.write(line)
    in_file.close()
    out_file.close()

    os.system('cp ' + code_dir + 'genTab/convertPairs.pl ' + data_dir + 'data/')
    os.system(data_dir + 'data/convertPairs.pl')
    os.system('rm ' + data_dir + 'data/convertPairs.pl')

    ###################################################
    ## Determine genTab parameters                   ##
    ###################################################
    max_neighbours = -10
    f = open(data_dir + 'data/protein_pairs_index.txt', 'r')
    for line in f:
        if len(line.strip().split()) - 2 > max_neighbours:
            max_neighbours = len(line.strip().split()) - 2
    f.close()

    ###################################################
    ## Set genTab parameters and recompile and run   ##
    ###################################################
    cmd = 'cp ' + code_dir + 'genTab/PIPE.h ' + code_dir + 'genTab/PIPE.h_bak'
    print cmd
    os.system(cmd)
    in_file  = open(code_dir + 'genTab/PIPE.h_bak', 'r')
    out_file = open(code_dir + 'genTab/PIPE.h', 'w')

    for line in in_file:
        if line.startswith('#define MAX_NUM_PROTEINS'):
            out_file.write('#define MAX_NUM_PROTEINS ' + str(len(proteins)) + '\n')
        elif line.startswith('#define MAX_PROTEIN_LEN'):
            # for some reason genTab throws an error unless we do + 1
            out_file.write('#define MAX_PROTEIN_LEN ' + str(longest_protein + 1) + '\n') 
        elif line.startswith('#define MAX_NEIGHBOURS'):
            out_file.write('#define MAX_NEIGHBOURS ' + str(max_neighbours) + '\n')
        else:
            out_file.write(line)
    in_file.close()
    out_file.close()

    # RECOMPILE
    print 'Compiling Gentab ...'
    cwd = os.getcwd()
    os.chdir(code_dir + 'genTab')
    os.system('make clean')
    os.system('make')

    os.chdir(cwd)    
    
    print 'Now run genTab, then next setup script, e.g:'
    print "srun -t 1:00:00 --ntasks 24 --mem-per-cpu=1G \\" 
    print "./code/genTab/genTab \\"
    print "./data/protein_sequences.txt \\"
    print "./database \\"
    print "./data/genTab_org.txt & disown \\"