#!/usr/bin/python 

import sys
import os
import glob
from random import shuffle

# Script accepts up to 2 arguments:
# Either:
#   ./setup.py          - to run both genTab and MP-PIPE
#   ./setup.py genTab   - to run genTab and gather data
#   ./setup.py MP-PIPE  - to run MP-PIPE and upload local database
#   ./setup.py test             -to run testing for both
#   ./setup.py genTab test      -run genTab and compare. Not expected and may produce error
#   ./setup.py MP-PIPE test     -run MP-PIPE and test output

if len(sys.argv) == 2 and sys.argv[1] == 'test':
    # run in testing mode
    testing = True
elif len(sys.argv) == 3 and sys.argv[2] == 'test':
    testing = True
else:
    testing = False

run_gentab = True
run_mppipe = True

if len(sys.argv) > 1 and sys.argv[1] == 'genTab':
    run_mppipe = False
if len(sys.argv) > 1 and sys.argv[1] == 'MP-PIPE':
    run_gentab = False


def get_num_hosts(hostfile):
    c = 0
    f = open(hostfile, 'r')
    for line in f:
        if not line.startswith('#'):
            c = c + 1
    f.close()
    return c


####################################################
## Check to see that these parameters are correct ##
####################################################
local_dir = '/home/bradbarnes/code/newSoyPIPE/' # this directory will contain the correct PIPE folder
remote_dir = '/home/bradbarnes/temp_pipe/'    # this is where you code will run out of remotely   

update = local_dir + 'PIPE/code/scripts/update.py'

collect = local_dir + 'PIPE/code/scripts/collect.py'

genTab = remote_dir + 'PIPE/code/genTab_New/genTab'
genTab_hostfile = local_dir + 'PIPE/code/genTab_New/genTab_hosts'
num_genTab_hosts = get_num_hosts(genTab_hostfile)

mp_pipe = remote_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2_New'
mp_pipe_hostfile = local_dir + 'PIPE/code/MP-PIPE2_New/PIPE_hosts'
num_mp_pipe_hosts = get_num_hosts(mp_pipe_hostfile)
print "%d pipe hosts" % num_mp_pipe_hosts

organism_name = 'SCN'


########################################################
## Check the sequence file to ensure it's formatted   ##
## correctly and that the proteins only contain valid ##
## amino acids (check_sequences_and_pairs.py)         ##
########################################################
cur_dir = os.getcwd()
os.chdir(local_dir + 'PIPE/data/' + organism_name)

print 'Checking protein sequences...'
AAs = ['A','R', 'N','D', 'C', 'Q','E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S','T', 'W', 'Y', 'V', 'B', 'Z', 'X', 'U']
proteins = set([])
longest_protein = 0;

sequence_file = open('data/protein_sequences.txt', 'r')
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

sequence_file.close()


#########################################################
## Check the known pairs file to ensure it's formatted ##
## correctly and that the pairs only contain proteins  ##
## from the previously processed sequence file         ##
## (check_sequences_and_pairs.py)                      ##
#########################################################
print 'Checking protein pairs...'
num_known_pairs = 0

pair_file = open('data/protein_pairs.txt', 'r')
for line in pair_file:
    cols = line.strip().split('\t')

    if len(cols) != 2:
        print 'ERROR: Protein pairs file not formatted properly'
        sys.exit(0)
    
    if cols[0] not in proteins:
        print 'ERROR: Protein ' + cols[0] + 'appears in protein pairs but not in sequence file'
        sys.exit(0)
    
    if cols[1] not in proteins:
        print 'ERROR: Protein ' + cols[1] + 'appears in protein pairs but not in sequence file'
        sys.exit(0)
    
    num_known_pairs = num_known_pairs + 1

pair_file.close()


if run_gentab:
    ################################################
    ## Create the interaction graph file from the ##
    ## protein pairs and sequence files           ##
    ################################################
    print 'Creating protein index file...'
    os.system('cp ' + local_dir + 'PIPE/code/genTab_New/convertPairs.pl ./')
    os.system('./convertPairs.pl')
    os.system('rm convertPairs.pl')


    ######################################################
    ## Create the directory to store the database files ##
    ######################################################
    os.system('mkdir database_New')


    #################################################
    ## Remake genTab is any changes have been made ##
    #################################################
    os.system('make -C ' + local_dir + 'PIPE/code/genTab_New/')

    #####################################################################
    ## Ensure all of the nodes on the cluster have the proper file     ##
    ## structure by rsyncing the current deiectory to all of the nodes ##
    #####################################################################
    os.system(update + ' ' + local_dir + 'PIPE ' + remote_dir + ' ' + genTab_hostfile)


    #####################################################################
    ## Run the program to generate all of the databwse files on all of ##
    ## the nodes. This step can take anywhere from 10 to 90 minutes    ##
    #####################################################################
    print ('mpirun -np ' + str(num_genTab_hosts*1-0) + ' -hostfile ' +  genTab_hostfile + ' ' +  genTab + ' ' +
            remote_dir + 'PIPE/data/' + organism_name +'/data/protein_sequences.txt ' +
            remote_dir + 'PIPE/data/' + organism_name +'/database_New ' + 
            remote_dir + 'PIPE/data/' + organism_name +'/data/genTab_org.txt '  + 
            remote_dir + 'PIPE/data/' + organism_name +'/data/protein_pairs_index.txt')
    os.system('mpirun -np ' + str(num_genTab_hosts*1-0) + ' -hostfile ' +  genTab_hostfile + ' ' +  genTab + ' ' + 
                remote_dir + 'PIPE/data/' + organism_name +'/data/protein_sequences.txt ' +
                remote_dir + 'PIPE/data/' + organism_name +'/database_New ' +
                remote_dir + 'PIPE/data/' + organism_name +'/data/genTab_org.txt ' +
                remote_dir + 'PIPE/data/' + organism_name +'/data/protein_pairs_index.txt')


    #####################################################################################
    ## Collect all of the database files that were created locally on all of the nodes ##
    #####################################################################################
    os.system(collect + ' ' + local_dir + ' ' + remote_dir + 'PIPE ' + genTab_hostfile)


###################################################
## Determine all other necessary PIPE parameters ##
###################################################
num_set_bits = 0
while num_set_bits < len(proteins):
    num_set_bits = num_set_bits + 64

max_neighbours = -10
f = open('data/protein_pairs_index.txt', 'r')
for line in f:
    if len(line.strip().split()) - 2 > max_neighbours:
        max_neighbours = len(line.strip().split()) - 2
f.close()

max_db_file = 0
db_files = glob.glob('./database_New/*')
for db_file in db_files:
    if os.path.getsize(db_file) > max_db_file:
        max_db_file = os.path.getsize(db_file)
#db_files = glob.glob('./database/*')
#for db_file in db_files:
#    if os.path.getsize(db_file) > max_db_file:
#        max_db_file = os.path.getsize(db_file)

#############################
## Write paramters to file ##
#############################
f = open('parameters.txt', 'w')

f.write('num proteins\t' + str(len(proteins)) + '\n')
f.write('num known pairs\t' + str(num_known_pairs) + '\n')
f.write('num pairs\t' + str((len(proteins) * (1+len(proteins)))/2) + '\n')
f.write('max_db_file\t' + str(max_db_file) + '\n')
f.write('num_set_bits\t' + str(num_set_bits) + '\n')
f.write('max_protein_len\t' + str(longest_protein) + '\n')
f.write('max_neighbours\t' + str(max_neighbours) + '\n')
f.close()

if run_mppipe:
    #############################################################################
    ## Update the mp-pipe2 code to reflect the new parameters and recompile it ##
    #############################################################################
    print 'cp ' + local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2.c ' + local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2.c_bak'
    os.system('cp ' + local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2.c ' + local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2.c_bak')
    in_file = open(local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2.c_bak', 'r')
    out_file = open(local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2.c', 'w')

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

    print ('mpicc -O3 -fopenmp -Wall ' + local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2.c -m64 -lm -o ' +
                        local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2_New')
    os.system('mpicc -O3 -fopenmp -Wall ' + local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2.c -m64 -lm -o ' + 
                local_dir + 'PIPE/code/MP-PIPE2_New/mp-pipe2_New')


    ###################################################
    ## Send the complete database and newly compiled ##
    ## mp-pipe back out to all of the nodes          ##
    ###################################################
    os.system(update + ' ' + local_dir + 'PIPE ' + remote_dir + ' ' + mp_pipe_hostfile)


    #############################################
    ## Create all-to-all input file            ##
    ## Format: first line = number of pairs    ##
    ## each other line: protein_ID\tprotein_ID ##
    #############################################
    os.system('mkdir input')
    os.system('mkdir output')

    if not testing:
        pass
        #os.system(local_dir + 'PIPE/code/scripts/create_all_to_all_input.py ' + str(len(proteins)) + 
        #    ' input/' + organism_name + '.in')
    else:
        # copy testing input file so that it is run
        #os.system('cp ' + local_dir + 'PIPE/data/' + organism_name + '/test/test.in input/' + organism_name + '.in')
        pass

    #################
    ## run MP-PIPE ##
    #################
    #mpirun -np <# processes> mp-pipe2 <protein pairs input file> <output file> <PIPE Interaction Graph File> <PIPE DB dir
    print ('mpirun -np ' + str(num_mp_pipe_hosts) + ' -hostfile ' + mp_pipe_hostfile + ' ' +
            mp_pipe + ' input/' + organism_name + '.in output/' + organism_name + '.out ' +
            remote_dir + 'PIPE/data/' + organism_name + '/data/protein_pairs_index.txt '+
            remote_dir + 'PIPE/data/' + organism_name + '/database_New/databaseSxH ' + 
            remote_dir + 'PIPE/data/' + organism_name + '/data/PIPE_org.txt')
    os.system('mpirun -np ' + str(num_mp_pipe_hosts) + ' -hostfile ' + mp_pipe_hostfile + ' ' + 
                mp_pipe + ' input/' + organism_name + '.in output/' + organism_name + '.out ' +
                remote_dir + 'PIPE/data/' + organism_name + '/data/protein_pairs_index.txt '+
                remote_dir + 'PIPE/data/' + organism_name + '/database_New/databaseSxH ' +
                remote_dir + 'PIPE/data/' + organism_name + '/data/PIPE_org.txt')


################################
## compare output if testing  ##
################################
if testing:
    os.system(local_dir + 'PIPE/code/scripts/check_output.py ' + 
        'output/' + organism_name + '.out ' +
        local_dir + 'PIPE/data/' + organism_name + '/test/test.out ')


os.chdir(cur_dir)            

