#!/usr/bin/python
#update local_dir remote_dir hostfile

import os
import sys

if len(sys.argv) != 4:
    print 'USAGE: update.py local_dir remote_dir hostfile'
    print '\nrsyncs local_dir to remote_dir on aall non-commented machines in hostfile'
    print 'WARNING: this uses --delete. Any data in remote_dir not in local_dir will be deleted\n'
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
    cmd = "ssh " + machine + " 'mkdir " + remote_dir + "'"
    print cmd
    os.system(cmd)

for machine in machines:
    #cmd = 'rsync -az --delete /aschoenr/ ' + machine + ':/aschoenr'
    cmd = ('rsync -az --delete --exclude=input --exclude=output --exclude=loocv --exclude=*.py ' + 
            '--exclude=create_roc --exclude=hello_world --exclude=scripts --exclude=landscapes ' + 
            local_dir + ' '  + machine + ':' + remote_dir)
    print cmd
    os.system(cmd)
