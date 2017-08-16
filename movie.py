#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import imp
from small_script.myFunctions import shrinkage

parser = argparse.ArgumentParser(description="written by Wei Lu.")
parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-p", "--plot", help="only plot mode",
                    action="store_true")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-l", "--last", action="store_true", default=False)
parser.add_argument("-m", "--mode",
                    type=int, default=2)
args = parser.parse_args()
# render Tachyon frame450.dat '/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86' -aasamples 12 %s -format TARGA -o frame450.tga -res 2000 2000

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

protein_name = args.protein.split('.')[0]

if args.mode == 1:
    # do("python3 ~/opt/small_script/delete_lammps_frame.py")
    shrinkage(shrink_size=1, max_frame=2000)
    # shrinkage()
    do("cp ../{}.seq .".format(protein_name))
    do("python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py small.lammpstrj movie "+protein_name+".seq")
    do("cp ~/opt/plot_scripts/2xov_movie_bicelle.tcl .")
    do("cp ~/opt/plot_scripts/movie_bicelle_no_smooth.tcl .")
    do("cp ~/opt/plot_scripts/movie.tcl .")
    if args.plot:
        do("/Applications/VMD\ 1.9.3.app/Contents/MacOS/startup.command -e 2xov_movie_bicelle.tcl")

if args.mode == 2:
    if not args.last:
        do("cp ../{}.seq .".format(protein_name))
        do("python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie "+protein_name+".seq")
        do("cp ~/opt/plot_scripts/2xov_movie_bicelle.tcl .")
        do("cp ~/opt/plot_scripts/movie_bicelle_no_smooth.tcl .")
        do("cp ~/opt/plot_scripts/2xov_movie_bicelle_no_smooth.tcl .")
        do("cp ~/opt/plot_scripts/movie.tcl .")
    else:
        x = os.listdir()
        y = [int(int(f.split(".")[1])/1000) for f in os.listdir() if f.startswith("restart.")]
        last_frame = max(y)
        print("last Frame: {}".format(last_frame))
        do("python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj last "+protein_name+".seq {}".format(last_frame))
        do("cp ~/opt/plot_scripts/last_frame.pml .")
        do("pymol last_frame.pml")
# if(args.number == -1):
#     # do("cp ~/opt/pulling/2xov.seq .")
#     if(not args.plot):
#         os.system("python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie "+protein_name+".seq")
#         # os.system(
#         #     "python2 ~/opt/script/BuildAllAtomsFromLammps.py \
#         #     dump.lammpstrj movie")
# else:
#     os.system("ghead -n 399648 dump.lammpstrj > part_dump")
#     if(not args.plot):
#         os.system("python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py dump.lammpstrj movie "+protein_name+".seq")
#         # os.system(
#         #     "python2 ~/opt/script/BuildAllAtomsFromLammps.py \
#         #     part_dump movie")
# os.system("cp ~/opt/plot_scripts/*.tcl .")
