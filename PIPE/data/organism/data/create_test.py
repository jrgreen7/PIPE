#!/usr/bin/python

import shutil
from random import randint

shutil.copyfile('protein_pairs.txt','test_input.names')

proteins = []
with open('protein_sequences.txt','r') as prots:
    for l in prots:
        proteins.append(l.strip().split('\t')[0])

pairs = set()
with open('protein_pairs.txt','r') as pairs_file:
    for l in pairs_file:
        c = l.strip().split('\t')
        k1 = c[0] + '\t' + c[1]
        k2 = c[1] + '\t' + c[0]
        pairs.add(k1)
        pairs.add(k2)
        
negSet = set()
with open('test_input.names','a+') as out:
    negs = 66084
    added = 0
    while added < negs:
        # get random int1 and 1
        p1 = randint(0,len(proteins)-1)
        p2 = randint(0,len(proteins)-1)
        # check not positive or already added
        key1 = proteins[p1] + '\t' + proteins[p2]
        key2 = proteins[p2] + '\t' + proteins[p1]
        if not (key1 in pairs or key2 in pairs or key1 in negSet or key2 in negSet):
            # add it and reverse to negSet and increment count
            negSet.add(key1)
            added += 1
            out.write(key1 + '\n')
            
