""" autosubmitter_cedar.py
Author: Kevin Dick
Date: 2020-05-08
---
Autosubmits jobs to the Cedar queue.
"""
import os
from time import sleep 

base_dir = '/home/williamm/scratch/Deep-PIPE-Sites/PIPE4/'
sub_dir  = base_dir + 'submissions/'
land_dir = base_dir + 'landscapes/'
SLEEP    = 5

for prot in [x.split('.')[0].strip() for x in os.listdir(sub_dir)]:
    land = land_dir + prot
    cmd  = 'sbatch ' + sub_dir + prot + '.sub'
    print('Running: ' + cmd + '\tin: ' + land)

    # Make the landscape directory & move into ensure PIPE-Sites landscapes end up here
    os.mkdir(land)
    os.chdir(land)

    # Enqueue the Job
    os.system(cmd)
    sleep(SLEEP)
