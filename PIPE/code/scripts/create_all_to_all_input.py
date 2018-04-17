#!/usr/bin/python
#create_all_to_all.py num_proteins out_filename

import sys
from random import shuffle

if len(sys.argv) != 3:
    print 'USAGE: create_all_to_all.py num_proteins out_filename'
    print '\ncreates an MP-PIPE input file for all possible ID pairs'
    print 'for num_proteins proteins. written to out_filename\n'
    sys.exit(0)
    
np = int(sys.argv[1])
fn = sys.argv[2]

pairs = []

f = open(fn, 'w')
f.write(str(int((np+1)*np/2.)) + '\n')

for i in range(np):
    for j in range(i, np):
        f.write(str(i) + '\t' + str(j) + '\n')

f.close()
