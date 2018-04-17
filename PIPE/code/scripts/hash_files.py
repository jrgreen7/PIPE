#!/usr/bin/python

import os
import sys

if len(sys.argv) != 4:
    print 'USAGE: hash_files.py remote_dir hostfile outputfil'
    sys.exit(0)

remote_dir = sys.argv[1]
hostfile = sys.argv[2]
output_file = sys.argv[3]

machines = []
f = open(hostfile, 'r')
for line in f:
    if not line.startswith('#'):
        machines.append(line.strip().split()[0])
f.close()

for machine in machines:
    cmd = 'ssh %s "sh -c \\\"(nohup find %s -type f -exec sha1sum "{}" + > %s) > /dev/null &\\\""' % (machine, remote_dir, output_file)
    print cmd
    os.system(cmd)

