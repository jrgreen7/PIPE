#!/usr/bin/python
#check_load.py hostfile

import os
import sys

hostfile = sys.argv[1]

machines = []
f = open(hostfile, 'r')
for line in f:
    if not line.startswith('#') and line.strip().split()[0] not in machines:
        machines.append(line.strip().split()[0])
f.close()

for machine in machines:
    print machine
    os.system('ssh '+ machine + ' uptime')
