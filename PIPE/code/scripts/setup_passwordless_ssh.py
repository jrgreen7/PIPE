#!/usr/bin/python

import sys
import os

os.system("cd ~")
os.system("ssh-keygen -t rsa")

for i in range(1, 6):
    cmd = ("cat .ssh/id_rsa.pub | ssh desk0" +
            str(i) + " 'mkdir .ssh && cat >> .ssh/authorized_keys'")
            #str(i) + " 'cat >> .ssh/authorized_keys'")
    os.system(cmd)

for i in range(1, 24):
    if i < 10:
        cmd = "cat .ssh/id_rsa.pub | ssh node0" 
    else:
        cmd = "cat .ssh/id_rsa.pub | ssh node"
    cmd = cmd + str(i) + " 'mkdir .ssh && cat >> .ssh/authorized_keys'"
    #cmd = cmd + str(i) + " 'cat >> .ssh/authorized_keys'"
    os.system(cmd)
