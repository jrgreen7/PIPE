#!/usr/bin/python
#collect.py local_dir remote_dir hostfile

import os
import sys

if len(sys.argv) != 4:
    print 'USAGE: collect.py local_dir remote_dir hostfile'
    print '\nrsyncs remote_dir on all non-commented machines in hostfile'
    print 'to local_dir. Used to collect data spread across the cluster\n'
    sys.exit(0)

local_dir = sys.argv[1]
remote_dir = sys.argv[2]
hostfile = sys.argv[3]

machines = []
f = open(hostfile)
for line in f:
    if not line.startswith('#') and line.strip().split()[0] not in machines:
        machines.append(line.strip().split()[0])
f.close()

for machine in machines:
    cmd = 'rsync -avz ' + machine + ':' + remote_dir + ' ' + local_dir
    print cmd
    os.system(cmd)
