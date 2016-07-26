#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp
import glob

# os.chdir('')
name_list = glob.glob('*.dat')
# os.chdir('..')
print(name_list)
# os.system("gnuplot -e 'out_file_name="test"; in_file_name="T120_400"'")
for name in name_list:
    print("gnuplot -e \"out_file_name='%s'; in_file_name='%s'\" \
        ~/opt/plot_scripts/list_q.gp" % (name[0:-4]+'.pdf', name))
    os.system("gnuplot -e \"out_file_name='%s'; in_file_name='%s'\" \
        ~/opt/plot_scripts/list_q.gp" % (name[0:-4]+'.pdf', name))
