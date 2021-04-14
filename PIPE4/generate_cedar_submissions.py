#!/usr/bin/python
""" generate_submissions.py
Author: Kevin Dick
Date: 2019-02-02
---
Creates individual job submission BASH scripts to 
submit to the biocluster.
These need to be SCP'ed to the biocluster for running.
"""
import os, sys

NUM_THREADS    = 10
MEM_PER_THREAD = '128g' #'130G'#'85G' #'153G' #'126G' #'87G' #'156G'
JOB_TIME       = '00:30:00'

job_id = 'cvd'#'ep4'#'hivhiv' #'hivhiv' #'hs-hiv' #'soyV4-human' #'soyV4-scnV2' #'soyV2-human' #'soyV2-soyV2' #'soyV2-scnV2'

data_dir = '/home/williamm/scratch/final-wrap-up/PIPE/PIPE4/'
code_dir = data_dir + 'code/'

# Read in the list of completed job
#completed = [f.replace('.out', '') for f in os.listdir(data_dir + 'output/')]

# Read in all of the input files
#inputs = [f.replace('.in', '') for f in os.listdir(data_dir + '3hr-input/')]
inputs = [f.replace('.in', '') for f in os.listdir(data_dir + 'input/')]

for i in inputs:
	# Skip those that are allready finished...
	#if i in completed: 
	#	print 'Already computed : ' + i
	#	continue

	print 'Creating Submission File: ' + i + '.sub'
	sub_f = data_dir + 'submissions/' + i + '.sub'

	# Start writing out the submission file for this input...
	f = open(sub_f, 'w')
	
	# Submission Header
	f.write('\n'.join(['#!/bin/bash',
                           '#SBATCH --account=def-jrgreen', # Ensure we use the Cedar Covid Allocation 
                           '#SBATCH --job-name=pipe4-' + i + '-o2a',
                           '#SBATCH --output=' + data_dir + 'logs/' + i + '.log',
                           '#SBATCH --mem=' + str(MEM_PER_THREAD),
                           '#SBATCH --nodes=1',
                           '#SBATCH --ntasks-per-node=' + str(NUM_THREADS),
                           '#SBATCH --cpus-per-task=1',
                           '#SBATCH --time=' + JOB_TIME,
                           '\n']))

	# The actual run command
	f.write(' '.join(['mpirun',
			  '-np ' + str(NUM_THREADS),
			  code_dir + 'MP-PIPE2/mp-pipe2',
                          data_dir + 'input/' + i + '.in',
			  data_dir + 'output/' + i + '.out',
			  data_dir + 'data/protein_pairs_index.txt',
			  data_dir + 'data/protein_pairs_indexes.txt',
			  data_dir + 'database',
			  data_dir + 'data/PIPE_org.txt\n']))
	f.close()
print 'Execution Complete!~'
