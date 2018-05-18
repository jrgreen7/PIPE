#!/usr/bin/python
#check_output.py PIPE_output_file1 PIPE_output_file2

import sys

if len(sys.argv) != 3:
    print 'USAGE: check_output.py PIPE_output_file1 PIPE_output_file2'
    print '\nloads both result files to memory, sets the runtime of each pair to 0'
    print 'and then ensures each record from file1 is in file 2 and vice versa\n'
    sys.exit(0)

f1 = open(sys.argv[1], 'r')
f2 = open(sys.argv[2], 'r')

a = set([])
b = set([])
errors = 0

for line in f1:
	cols = line.strip().split('\t')
	cols[4] = 0
	a.add(str(cols))

for line in f2:
	cols = line.strip().split('\t')
	cols[4] = 0
	b.add(str(cols))
	
f1.close
f2.close

print "CHECKING FILE1 AGAINST FILE2"

for line_a in a:
	if line_a not in b:
		print line_a
		errors = errors + 1
	
print str(errors) + " errors."

errors = 0

#print "CHECKING FILE2 AGAINST FILE1"
#
#for line_b in b:
#	if line_b not in a:
#		print line_b
#		errors = errors + 1
#	
#print str(errors) + " errors."
#
