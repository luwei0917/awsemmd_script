#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import subprocess
# compute cross Q for every pdb pair in one folder
cross_q = open('cross_q', 'w')

n = 20
for i in range(1, n+1):
    for j in range(1, n+1):
        script = "python2 ~/opt/small_script/CalcQValueFrom2Pdb.py {}.pdb {}.pdb".format(i, j)
        # x = subprocess.check_output()
        result = subprocess.check_output(script, shell=True).decode("utf-8")
        # print(result.strip('\n'))
        cross_q.write(result.strip('\n'))
        cross_q.write(" ")
    cross_q.write("\n")
