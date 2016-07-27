#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp
import glob
parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        generate the fragment memory from lowest energy frame")

memfile = open("lowest_energy.mem", 'w')
memfile.write("[Target]" + "\n")
memfile.write("query" + "\n" + "\n")

memfile.write("[Memories]" + "\n")

os.chdir("lowest_energy")
filelist = glob.glob("*.pdb")


for frag in filelist:
    # print(frag.split(".")[0])
    name = frag.split(".")[0]
    os.system("python2 ~/opt/script/Pdb2Gro.py %s.pdb %s.gro" % (name, name))
    #./fraglib/4a2ba.gro 1 86 9 1
    length = int((sum(1 for line in open(name+".gro"))-2)/5)
    # print(length)
    memfile.write("./lowest_energy/" + name + ".gro 2 2 "+str(length)+" 1\n")
os.chdir("..")
memfile.close()
