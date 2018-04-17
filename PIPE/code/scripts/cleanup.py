#!/usr/bin/python
#cleanup.py remote_dir hostfile

import os
import sys

if len(sys.argv) != 3:
    print 'USAGE: clearnup.py remote_dir hostfile'
    print '\nremoves (rm -r) remote_dir on all non-commented machines in hostfile\n'
    sys.exit(0)

remote_dir = sys.argv[1]
hostfile = sys.argv[2]

machines = []
f = open(hostfile, 'r')
for line in f:
    if not line.startswith('#'):
        machines.append(line.strip().split()[0])
f.close()

for machine in machines:
    cmd = "ssh " + machine + " 'rm -r " + remote_dir + "'"
    print cmd
    os.system(cmd)
