# PIPE4 Documentation for Compute Canada

Date: 2020-10-06
Author: Kevin Dick
Environment: Cedar (Compute Canada)

---

The Protein Interaction Prediction Engine (PIPE) was originally developed in 2006 and refined through various versions over the coourse of a decade. In each iteration, the algorithm was optimized for time (computed pair/second), space (amount of RAM), and parrallelization (high performance computing implementation as "embarrassingly parallel").

The latest iteration of the PIPE algorithm, PIPE4, is adapted to flexibly support inter- and cross-species prediction schemas. Consequently, the PIPE4 algorithm requires special consideration in how the configuration files are setup.

## PIPE4 File Structure

The typical PIPE4 file structure on a Compute Canada cluster such as Cedar will have the following subdirectories:
```
code\
data\
database\  
input\
output\
landscapes\  
submissions\
logs\
```

The `code\` directory contains two subdirectories, the `genTab` directory with the precomputation codebase and the `MP-PIPE2` directory with the latest PIPE4 codebase.

The `data\` directory must contain the `protein_pairs.txt`, `protein_sequences.txt`, `genTab_org.txt`, and `PIPE_org.txt` files, prepared according to the specifications.

The `database\` directory will contain two database files for eech of the precomputed proteins according to the `genTab_org.txt` file. If genTab fails, these files will be corrupt and you will need to restart the entire job.

The `input\` directory will contain all of the "input" files to determine what predictions need to be made. Note: there is a difference in `protein_sequences.txt` numbering (1-indexed) and `protein_sequence_index(s).txt` (0-indexed); be very careful not to have 'off-by-one' errors here. These files need to be generated according to the protein indecies.

The `output\` directory is where all the PIPE4 output files will be written.

The `landscapes\` directory s where the landscapes are written.

Compute Canada requires job submission using submission scripts. The `submissions\` directory is meant to help organize the individual scripts to submit jobs to the Compute Canada queue.

The `logs\` directory is where the inidividual log files for jobs may be written for monitoring.

