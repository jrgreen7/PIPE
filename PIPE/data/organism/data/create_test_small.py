#!/usr/bin/python

import shutil
from random import randint


num_prots = 1000
num_test = 10000

with open('protein_sequences.txt_full','r') as inf, open('protein_sequences.txt','w') as out:
    for i in range(num_prots):
        out.write(inf.readline())

proteins = []
with open('protein_sequences.txt','r') as prots:
    for l in prots:
        proteins.append(l.strip().split('\t')[0])

proteins_set = set(proteins)
pairs = set()
with open('protein_pairs.txt_full','r') as inf, open('protein_pairs.txt','w') as out:
    for l in inf:
        c = l.strip().split('\t')
        if c[0] in proteins_set and c[1] in proteins_set:
            out.write(l)
            k1 = c[0] + '\t' + c[1]
            k2 = c[1] + '\t' + c[0]
            pairs.add(k1)
            pairs.add(k2)
        

shutil.copyfile('protein_pairs.txt','test_input.names')
negSet = set()
with open('test_input.names','a+') as out:
    added = 0
    while added < num_test:
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
            
