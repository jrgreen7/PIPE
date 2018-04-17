#!/usr/bin/python
#python set_PIPE_cutoff.py input_file output_file cutoff

import sys

if len(sys.argv) != 4:
    print 'USAGE: set_PIPE_cutoff.py input_file output_file cutoff'
    print '\nreads the input_file (produced by PIPE) and writes all'
    print 'records with a SW score greater than or equal to cutoff\n'
    sys.exit(0)

f1 = open(sys.argv[1], 'r')
f2 = open(sys.argv[2], 'w')
t = float(sys.argv[3])

for line in f1:
    cols = line.split('\t')
    if(float(cols[5]) >= t):
        f2.write(line)
	
f1.close()
f2.close()
	
