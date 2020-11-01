""" create_input_files_exclusions.py
Authro: Kevin Dick
Date: 2018-12-02
---
Creates each of the input files for the AGCAN runs.
All input files are produced with respect to one reference organism.
NOTE: LOADS A STATUS.DB FILE FOR EACH JOB WITH THE SPECIFIC EXCLUSIONS...
"""
from multiprocessing import Pool, cpu_count

org_dir = '/home/williamm/scratch/Deep-PIPE-Sites/PIPE4/'
SKIP_BOOL = False
SKIP_INDICIES = [3902, 4675, 4676, 4757]

def create_input_files(job):
        print('Preparing input files for job ' + str(job[0]))
        job_id, orgA_start, orgA_end, orgB_start, orgB_end  = job
    
        # Setup the paths to files and dirs...
        input_dir = str(org_dir) + '/input/'
        prot_seqs = str(org_dir) + '/data/protein_sequences.txt'

        # Load all lines in protein sequence to get protein name for inputs.
        print('Reading in protein sequences file...')
        f = open(prot_seqs, 'r')
        lines = f.readlines()
        f.close()

        # Iterate over the range of indecies for organismA and write out the protein input file
        for i in range(orgA_start - 1, orgA_end): # -1 for zero-indexing
                protein_name = lines[i].split('\t')[0].strip()
                output_file  = input_dir + protein_name + '.in'
                
                f = open(output_file, 'w')
                # First write out the number of predictions: (orgB_end - orgB_start)
                if SKIP_BOOL: f.write(str(orgB_end - orgB_start + 1 - len(SKIP_INDICIES)) + '\n')
                else: f.write(str(orgB_end - orgB_start + 1) + '\n') # +1 for the zero-indexing offset
                for j in range(orgB_start - 1, orgB_end): # -1 for zero-indexing
                    if SKIP_BOOL and  j in SKIP_INDICIES: continue
                    f.write(str(i) + '\t' + str(j) + '\n')
                f.close()          
        print('Completed for ' + str(job_id) + '!')

# RANGES ARE UNCHANGED...
#         <job_id>  ,<orgA_start>, <orgA_end>, <orgB_start>, <orgB_end>
#jobs = [('hs-hiv', 1, 23594, 23595, 23603)]
#jobs = [('hiv-hiv', 23595, 23603, 23595, 23603)]
#jobs  = [('hs-hs', 1, 23594, 1, 23594)]
#jobs = [('ep4', 1, 23594, 23595, 23606)]
#jobs =  [('cvd', 689, 705, 709, 21071)]
jobs =  [('cvd', 35, 705, 706, 21071)]

# PARALLEL
#with  Pool(cpu_count() - 2) as p: p.map(create_input_files, jobs)

# SERIAL
for j in jobs: create_input_files(j)
print("Execution Complete")

