#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp
import glob

os.chdir('data')
name_list = glob.glob('*')
os.chdir('..')
print(name_list)
# os.system("gnuplot -e 'out_file_name="test"; in_file_name="T120_400"'")
for name in name_list:
    print("gnuplot -e \"out_file_name='%s'; in_file_name='%s'\" \
        ~/opt/plot_scripts/list_q.gp" % ('plot/'+name+'.pdf', 'data/'+name))
    os.system("gnuplot -e \"out_file_name='%s'; in_file_name='%s'\" \
        ~/opt/plot_scripts/list_q.gp" % ('plot/'+name+'.pdf', 'data/'+name))
