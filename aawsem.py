#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
folder_list = ["T0792", "T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
# folder_list = ["T0815", "T0778"]
os.system("mkdir -p aawsemDec25")
os.chdir("aawsemDec25")
for protein_name in folder_list:
    os.system("mkdir -p "+protein_name)
    os.chdir(protein_name)
    os.system("mkdir -p "+protein_name)
    os.chdir(protein_name)
    os.system("cp ~/Research/project_AAWSEM/aasem/_output/"+protein_name+".fasta .")
    os.system("create_project.py " + protein_name)
    os.chdir("../..")
# os.system("stride %s.pdb > ssweight.stride" % protein_name)
# os.system("python2 ~/opt/script/stride2ssweight.py > ssweight")
# os.system("python2 ~/opt/script/GetCACADistancesFile.py %s native.dat" % protein_name)
# os.system("python2 ~/opt/script/GetCACoordinatesFromPDB.py %s nativecoords.dat" % protein_name)
# os.system("cp ~/opt/database/* .")
# os.system("python2 ~/opt/script/prepFragsLAMW_index.py \
#     cullpdb_pc80_res3.0_R1.0_d160504_chains29712 %s.fasta 20 0" % protein_name)
# os.system("cp ~/opt/AASEM/parameter/* .")
# os.system("cp ~/opt/parameter/* .")
# os.system("python2 ~/opt/script/Pdb2Gro.py %s.pdb %s.gro" % (protein_name, protein_name))
# os.system("cp ~/opt/variables.dat .")
# print("You need zim")
# grep -E "CB|CA  GLY" 1qjp.pdb > cbs.data
#awk '{if($8>12.7) print "1"; else if($8<-12.7) print "3";else print "2" }' cbs.data  > zim
#awk '{if($8>14.4) print "1"; else if($8<-14.4) print "3";else print "2" }'  temp.data  > zim
# os.system("echo '%s A' > pdbidlist" % protein_name)
# os.system("getTransmembraneRegions.py pdbidlist")
# os.system("mv %sA.tm zim" % protein_name)


# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

# parser = argparse.ArgumentParser(
#         description="This is a python3 script to\
#         automatically analysis the simulation")
#
# parser.add_argument("template", help="the name of template file")
# parser.add_argument("-n", "--number", type=int, default=20,
#                     help="Number of simulation run")
# parser.add_argument("-m", "--movie", type=int, default=-1,
#                     help="generate the movie,defalut is all")
# parser.add_argument("-p", "--plotOnly", help="only generate the plot",
#                     action="store_true")
# parser.add_argument("-s", "--steps", type=int, default=4,
#                     help="Simulation steps in unit of million,\
#                     default is 4 million, -1 means test run")
# args = parser.parse_args()
# list_of_max_q = []
#
# n = args.number
# if args.steps == -1:
#     n = 1  # also set n to be 1
#
# protein_name = args.template.strip('/')

# os.system("mkdir -p results")

# os.system("")
# for i in range(n):
#     # analysis
#     os.system("mkdir -p analysis/"+str(i))
#     os.chdir("analysis/"+str(i))
#     if not args.plotOnly:
#         # move necessary file into analysis folder
#         sys.stdout = open("chosen.txt", "w")
#         os.system("mv ../../simulation/"+str(i)+"/dump.lammpstrj .")
#         os.system("mv ../../simulation/"+str(i)+"/wham.dat .")
#         os.system("mv ../../simulation/"+str(i)+"/energy.dat .")
#         record_time = 0
#         with open('wham.dat') as input_data:
#             # Skips text before the beginning of the interesting block:
#             record_time = 0
#             max_q = 0
#             last_q = 0
#             next(input_data)
#             for line in input_data:
#                 # time, q = line.strip().split()
#
#                 time = int(line.strip().split()[0])
#                 q = float(line.strip().split()[1])
#                 if(q > max_q):
#                     record_time = time
#                     max_q = q
#                 last_q = q
#             list_of_max_q += [(max_q, record_time, last_q)]
#         time_step = record_time
#
#         print('ITEM: TIMESTEP')
#         with open('dump.lammpstrj') as input_data:
#             # Skips text before the beginning of the interesting block:
#             for line in input_data:
#                 if line.strip() == str(time_step):
#                     print(line.strip())  # Or whatever test is needed
#                     break
#             # Reads text until the end of the block:
#             for line in input_data:  # This keeps reading the file
#                 if line.strip() == 'ITEM: TIMESTEP':
#                     break
#                 print(line.strip())
#         sys.stdout.close()
#     if(args.movie == -1 or args.movie == i):
#         os.system(
#             "python2 ~/opt/script/BuildAllAtomsFromLammps.py \
#             dump.lammpstrj movie")
#
#     sys.stdout = open("final.txt", "w")
#     print('ITEM: TIMESTEP')
#     time_step = args.steps*1000*1000
#     with open('dump.lammpstrj') as input_data:
#         # Skips text before the beginning of the interesting block:
#         for line in input_data:
#             if line.strip() == str(time_step):
#                 print(line.strip())  # Or whatever test is needed
#                 break
#         # Reads text until the end of the block:
#         for line in input_data:  # This keeps reading the file
#             if line.strip() == 'ITEM: TIMESTEP':
#                 break
#             print(line.strip())
#     sys.stdout.close()
#
#     os.system(
#         "python2 ~/opt/script/BuildAllAtomsFromLammps.py \
#         chosen.txt chosen")
#     # os.system("cp ~/opt/plot_scripts/energy.plt .")
#     os.system(
#         "python2 ~/opt/script/BuildAllAtomsFromLammps_seq.py \
#         final.txt final.pdb ../../" +
#         protein_name+"/"+protein_name+".seq "+str(args.steps*1000))
#     # plots
#     os.system("cp ~/opt/plot_scripts/*.plt .")
#     os.system("cp ~/opt/plot_scripts/*.pml .")
#     os.system("/usr/local/bin/pymol -qc -r print_final.pml")
#     os.system(  # replace PROTEIN with pdb name
#             "sed -i.bak 's/PROTEIN/'" +
#             protein_name +
#             "'/g' show_origin.pml")
#
#     os.system(  # replace PROTEIN with pdb name
#             "sed -i.bak 's/NUMBER/'" +
#             str(i) +
#             "'/g' energy.plt")
#     os.system("gnuplot energy.plt")
#
#     os.system(  # replace PROTEIN with pdb name
#             "sed -i.bak 's/NUMBER/'" +
#             str(i) +
#             "'/g' q_value.plt")
#     os.system("gnuplot q_value.plt")
#
#     os.system("gnuplot detail_energy.plt")
#     # subprocess.Popen("gnuplot q_value.plt", env=my_env)
#     os.system("cp ~/opt/plot_scripts/*.tcl .")
#     os.system(  # replace PROTEIN with pdb name
#             "sed -i.bak 's/PROTEIN/'" +
#             protein_name +
#             "'/g' membrane_show.tcl")
#     os.system(  # replace PROTEIN with pdb name
#             "sed -i.bak 's/PROTEIN/'" +
#             protein_name +
#             "'/g' show.tcl")
#     os.system("cp ../../"+protein_name+"/*.pdb .")
#     os.system(
#             "python2 ~/opt/script/CalcRMSD.py "+protein_name+" \
#             dump.lammpstrj rmsd")
#     os.chdir("../..")
# if not args.plotOnly:
#     sys.stdout = open("analysis/list_of_max_q", "w")
#     for q in list_of_max_q:
#         print(q[0], q[1], q[2])  # max q, timestep of max q, last q
#     sys.stdout.close()
