#!/usr/bin/python 

import sys
import os

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
local_dir = '/storage/aschoenr/' # this directory will contain the correct PIPE folder
remote_dir = '/tmp/aschoenr/'    # this is where you code will run out of remotely   

update = local_dir + 'PIPE/code/scripts/update.py'

mp_pipe = remote_dir + 'PIPE/code/MP-PIPE2/mp-pipe2_loocv'
mp_pipe_hostfile = local_dir + 'PIPE/code/MP-PIPE2/PIPE_hosts'
num_mp_pipe_hosts = get_num_hosts(mp_pipe_hostfile)

organism_name = 's.cerevisiae'


#######################################
## create the loocv folder structure ##
#######################################
print 'creating loocv folder structure...'
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/input')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/output')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/cutoffs')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/roc_curves')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/stripped_output')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/cutoffs')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/roc_curves')
os.system('mkdir ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/stripped_output')


############################################
## create PIPE input file for known pairs ##
############################################
print 'creating positive and negative input files...'
os.system(local_dir + 'PIPE/code/scripts/convert_names_to_IDs.py ' + 
           local_dir + 'PIPE/data/' + organism_name + '/data/protein_pairs.txt ' +
           local_dir + 'PIPE/data/' + organism_name + '/data/protein_pairs_index.txt temp')

pairs = []
f = open('temp', 'r')
for line in f:
    pairs.append(line)
f.close()

f = open(local_dir + 'PIPE/data/' + organism_name + '/loocv/input/p.in', 'w')
f.write(str(len(pairs)) + '\n')
for pair in pairs:
   f.write(pair)
f.close()

os.system('rm temp')


#################################################
## create PIPE input files for  negative pairs ##
#################################################
num_proteins = 0
f = open(local_dir + 'PIPE/data/' + organism_name + '/data/protein_sequences.txt', 'r')
for line in f:
    num_proteins = num_proteins + 1
f.close()

for i in range(1,4):
    os.system(local_dir + 'PIPE/code/create_roc/create_random_ID_pairs.py ' +
                local_dir + 'PIPE/data/' + organism_name + '/loocv/input/p.in ' + 
                str(num_proteins) +  ' 100000 ' + 
                local_dir + 'PIPE/data/' + organism_name + '/loocv/input/n' + str(i) + '.in')


#############################################################################
## Update the mp-pipe2 code to reflect the new parameters and recompile it ##
#############################################################################
print 'updating and recompiling MP-PIPE code...'
in_file = open(local_dir + 'PIPE/code/MP-PIPE2/mp-pipe2.c', 'r')
out_file = open(local_dir + 'PIPE/code/MP-PIPE2/mp-pipe2_loocv.c', 'w')

for line in in_file:
    if line.startswith('//#define ENABLE_LOOCV'):
        out_file.write('#define ENABLE_LOOCV\n')
    elif line.startswith('static const int PACKET_SIZE = '):
        out_file.write('static const int PACKET_SIZE = 2500;\n')
    else:
        out_file.write(line)
in_file.close()
out_file.close()

print ('mpicc -O3 -fopenmp -Wall ' + local_dir + 'PIPE/code/MP-PIPE2/mp-pipe2_loocv.c -m64 -lm -o ' +
                    local_dir + 'PIPE/code/MP-PIPE2/mp-pipe2_loocv')
os.system('mpicc -O3 -fopenmp -Wall ' + local_dir + 'PIPE/code/MP-PIPE2/mp-pipe2_loocv.c -m64 -lm -o ' +
            local_dir + 'PIPE/code/MP-PIPE2/mp-pipe2_loocv')


#################################################################
## send the newly compiled code out to the rest of the cluster ##
## and run it on the positive and negative pairs               ##
#################################################################
os.system(update + ' ' + local_dir + 'PIPE ' + remote_dir + ' ' + mp_pipe_hostfile)

os.system('mpirun -np ' + str(num_mp_pipe_hosts) + ' -hostfile ' + mp_pipe_hostfile + ' ' + 
            mp_pipe + ' ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/input/p.in ' +
            local_dir + 'PIPE/data/' + organism_name + '/loocv/output/p.out '  +
            remote_dir + 'PIPE/data/' + organism_name + '/data/protein_pairs_index.txt '+
            remote_dir + 'PIPE/data/' + organism_name + '/database')

for i in range(1,4):
    os.system('mpirun -np ' + str(num_mp_pipe_hosts) + ' -hostfile ' + mp_pipe_hostfile + ' ' +
                mp_pipe + ' ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/input/n' + str(i) + '.in ' + 
                local_dir + 'PIPE/data/' + organism_name + '/loocv/output/n' + str(i) + '.out '  +
                remote_dir + 'PIPE/data/' + organism_name + '/data/protein_pairs_index.txt '+
                remote_dir + 'PIPE/data/' + organism_name + '/database')


#################################################
## process output files for further processing ##
#################################################
print 'processing PIPE output files...'
os.system('cut -f1,2,3 ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/output/p.out > ' +
            local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/stripped_output/p.out')
os.system('cut -f1,2,6 ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/output/p.out > ' + 
            local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/stripped_output/p.out')
for i in range(1,4):
    os.system('cut -f1,2,3 ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/output/n' + str(i) + '.out > ' + 
                local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/stripped_output/n' + str(i) + '.out')
    os.system('cut -f1,2,6 ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/output/n' + str(i) + '.out > ' + 
                local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/stripped_output/n' + str(i) + '.out')


#############################################################
## copy the 'roc' program and run it on all negative files ##
#############################################################
print 'creating roc curve files...'
os.system('cp ' + local_dir + 'PIPE/code/create_roc/roc ./')
#Usage: ./roc <pos file> <#pos> <neg file> <#neg>

for i in range(1,4):
    os.system('./roc ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/stripped_output/p.out ' + 
                str(len(pairs)) + ' ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/stripped_output/n' +
                str(i) + '.out 100000  > ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/cutoffs/n' + 
                str(i) + '.cutoffs')
    os.system('mv results.roc ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/roc_curves/n' + str(i) + '.roc')

    os.system('./roc ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/stripped_output/p.out ' + 
                str(len(pairs)) + ' ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/stripped_output/n' +
                str(i) + '.out 100000  > ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/cutoffs/n' + 
                str(i) + '.cutoffs')
    os.system('mv results.roc ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/roc_curves/n' + str(i) + '.roc') 
    
os.system('rm ./roc')


##################################################################
## copy the gunplot script to CWD, run it to produce roc curves ##
##################################################################
print 'creating roc curve...'

cur_dir = os.getcwd()
os.chdir(local_dir + 'PIPE/data/' + organism_name + '/loocv/')
os.system('cp ' + local_dir + 'PIPE/code/create_roc/loocv.plt ./')
os.system('gnuplot ./loocv.plt')
#os.system('mv loocv.eps ' + local_dir + 'PIPE/data/' + organism_name + '/loocv/')
#os.system('rm loocv.plt')
os.chdir(cur_dir)


#################################################################
## create cutoffs files for easy reference (99.95% specificity ##
#################################################################
print 'creating cutoff files...'
cutoff_file = open(local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/cutoffs/cutoffs.txt', 'w')
for i in range(1,4):
    f = open(local_dir + 'PIPE/data/' + organism_name + '/loocv/PIPE_score/cutoffs/n' + str(i) + '.cutoffs', 'r')
    line = f.readline()
    while float(line.split(',')[1].split(': ')[1]) < 0.99945:
        line = f.readline()
    f.close()
    cutoff_file.write('n' + str(i) + ': ' + line)
cutoff_file.close()

cutoff_file = open(local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/cutoffs/cutoffs.txt', 'w')
for i in range(1,4):
    f = open(local_dir + 'PIPE/data/' + organism_name + '/loocv/SW_score/cutoffs/n' + str(i) + '.cutoffs', 'r')
    line = f.readline()
    while float(line.split(',')[1].split(': ')[1]) < 0.99945:
        line = f.readline()
    f.close()
    cutoff_file.write('n' + str(i) + ': ' + line)
cutoff_file.close()


