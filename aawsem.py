#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import myPersonalFunctions

parser = argparse.ArgumentParser(description="This is my aawsem project headquarter")
parser.add_argument("--dec25", help="Run code on Dec 25", action="store_true", default=False)
parser.add_argument("--jan03", help="Run code on Jan 03 ", action="store_true", default=False)
parser.add_argument("--jan07", help="Run code on Jan 07 ", action="store_true", default=False)
parser.add_argument("--jan10", help="Run code on Jan 10 ", action="store_true", default=False)
parser.add_argument("--jan11", help="Run code on Jan 11 ", action="store_true", default=False)
parser.add_argument("--jan12", action="store_true", default=False)
parser.add_argument("--jan15", action="store_true", default=False)
parser.add_argument("-t", "--test", help="Test run", action="store_true", default=False)
parser.add_argument("--fix", action="store_true", default=False)
parser.add_argument("--move", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--may04", action="store_true", default=False)
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

if(args.may04):
    do("mkdir tertiary_T0782")
    do("cp ~/opt/AAWSEM/T0782.fasta .")
if(args.move):
    folder_list = ["T0766", "T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0766"]
    cd("aawsemJan15")
    for protein_name in folder_list:
        cd(protein_name)
        do("lowest_energy.py {}".format(protein_name))
        do("python2 ~/opt/LammpsPDBToCoordinates.py global_lowest_energy v1")
        do("python2 ~/opt/script/CoordinatesToWorkLammpsDataFile.py v1.coord data.v1")
        # do("cp -r {0} {0}_v1".format(protein_name))
        # do("mv simulation simulation_iteration_1")
        # do("~/opt/script/PdbCoords2Lammps.sh global_lowest_energy v1")
        do("cp data.v1 {0}/data.{0}".format(protein_name))
        do("mv list_of_lowest_potential_energy list_of_lowest_potential_energy_v1")
        do("run.py -n 20 -o "+protein_name+"\/")
        cd("..")
    cd("..")

if(args.fix):
    folder_list = ["T0766", "T0792", "T0778", "T0782", "T0833", "T0844"]
    os.chdir("aawsemJan15")
    for protein_name in folder_list:
        os.chdir(protein_name)

        os.system("cp {0}/frag1.mem {0}/frag.mem".format(protein_name))
        os.system("rm -r simulation")
        os.system("run.py -n 20 -o "+protein_name+"\/")
        os.chdir("..")

if(args.test):
    folder_list = ["T0792", "T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    for protein_name in folder_list:
        cd(protein_name)
        protein_size = myPersonalFunctions.file_width(protein_name+".seq")
        print(protein_size)
        n = 0
        nn = 0
        with open('Hybrid.mem', 'w') as w:
            with open('HE.mem', 'r') as f:
                he_lines = f.readlines()
            # print(ha_lines)
            with open('HO.mem', 'r') as f:
                ho_lines = f.readlines()
            for i in range(1, protein_size-7):
                fragGroup = []
                for line in ho_lines:
                    # _, loc, _, fragLens, _ = line.split(" ")
                    window, _, score, evalue, name, loc, j_start, fragLens, *weight = line.split()
                    locEnd = int(loc) + int(fragLens)
                    # print(loc, fragLens, locEnd)
                    if(i <= int(loc) and i+9 >= locEnd):
                        fragGroup += [line]
                groupLens = len(fragGroup)
                if(groupLens != 0):
                    nn += 1
                    for j in range(20):
                        w.write(fragGroup[j % groupLens])
                else:
                    n += 1
                    for j in range(20):
                        try:
                            w.write(he_lines[(i-1)*20+j])
                        except:
                            print(i*20+j)
        cd("..")
    # folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    #
    # print(folder_list)
    # # os.chdir("aawsemJan12")
    # for protein_name in folder_list:
    #     os.chdir(protein_name)
    #     os.chdir(protein_name)
    #     protein_size = myPersonalFunctions.file_len("ssweight")
    #     # print(protein_name, protein_size)
    #     with open("info.dat", 'w') as info:
    #         n = 0
    #         nn = 0
    #         with open('frag1.mem', 'w') as w:
    #             with open('ha_frag.mem', 'r') as f:
    #                 for i in range(4):
    #                     w.write(next(f))
    #                 ha_lines = f.readlines()
    #             # print(ha_lines)
    #             with open('ho_frag.mem', 'r') as f:
    #                 for i in range(4):
    #                     next(f)
    #                 ho_lines = f.readlines()
    #             for i in range(1, protein_size-7):
    #                 fragGroup = []
    #                 for line in ho_lines:
    #                     _, loc, _, fragLens, _ = line.split(" ")
    #                     locEnd = int(loc) + int(fragLens)
    #                     # print(loc, fragLens, locEnd)
    #                     if(i <= int(loc) and i+9 >= locEnd):
    #                         fragGroup += [line]
    #                 groupLens = len(fragGroup)
    #                 if(groupLens != 0):
    #                     nn += 1
    #                     for j in range(20):
    #                         w.write(fragGroup[j % groupLens])
    #                 else:
    #                     n += 1
    #                     for j in range(20):
    #                         try:
    #                             w.write(ha_lines[(i-1)*20+j])
    #                         except:
    #                             print(i*20+j)
    #             info.write("Use "+str(nn)+" HO and "+str(n)+" HA\n")
    #     os.chdir("../..")
    # os.chdir("aawsemJan12")
    # folder_list1 = ["T0792", "T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    # folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # # folder_list = ["T0792"]
    # # folder_list = list(set(folder_list1) - set(folder_list2))
    # for protein_name in folder_list:
    #     os.system("echo '"+protein_name+"' >> info.dat")
    #     os.system("cat {0}/{0}/info.dat >> info.dat".format(protein_name))

if(args.jan15):
    folder_list = ["T0766", "T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    os.system("mkdir -p aawsemJan15")
    os.chdir("aawsemJan15")
    for protein_name in folder_list:
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        os.system("cp -r ../../aawsemDec25/{0}/{0} .".format(protein_name))
        os.chdir(protein_name)
        os.system("cp /scratch/wl45/HE_Mem/{}/frag.mem he_frag_old.mem".format(protein_name))
        os.system("sed 's/\/work\/pw8\/mc70\/script/\/home\/wl45/g' he_frag_old.mem > he_frag.mem")
        os.system("cp frag.mem ho_frag.mem")
        protein_size = myPersonalFunctions.file_len("ssweight")
        # print(protein_name, protein_size)
        with open("info.dat", 'w') as info:
            n = 0
            nn = 0
            with open('frag1.mem', 'w') as w:
                with open('he_frag.mem', 'r') as f:
                    for i in range(4):
                        w.write(next(f))
                    ha_lines = f.readlines()
                # print(ha_lines)
                with open('ho_frag.mem', 'r') as f:
                    for i in range(4):
                        next(f)
                    ho_lines = f.readlines()
                for i in range(1, protein_size-7):
                    fragGroup = []
                    for line in ho_lines:
                        _, loc, _, fragLens, _ = line.split(" ")
                        locEnd = int(loc) + int(fragLens)
                        # print(loc, fragLens, locEnd)
                        if(i <= int(loc) and i+9 >= locEnd):
                            fragGroup += [line]
                    groupLens = len(fragGroup)
                    if(groupLens != 0):
                        nn += 1
                        for j in range(20):
                            w.write(fragGroup[j % groupLens])
                    else:
                        n += 1
                        for j in range(20):
                            try:
                                w.write(ha_lines[(i-1)*20+j])
                            except:
                                print(i*20+j)
                info.write("Use "+str(nn)+" HO and "+str(n)+" HE\n")
                print("Use "+str(nn)+" HO and "+str(n)+" HE\n")
        os.chdir("../..")

if(args.jan12):
    folder_list1 = ["T0792", "T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    folder_list2 = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    folder_list = list(set(folder_list1) - set(folder_list2))
    print(folder_list)
    # exit()
    os.system("mkdir -p aawsemJan12")
    os.chdir("aawsemJan12")
    for protein_name in folder_list:
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        os.system("cp -r ../../aawsemDec25/{0}/{0} .".format(protein_name))
        os.chdir(protein_name)
        os.system("cp ~/Research/CASP11/HA_Mem/{}/frag.mem ha_frag_old.mem".format(protein_name))
        os.system("sed 's/\/work\/pw8\/mc70\/script/\/home\/wl45/g' ha_frag_old.mem > ha_frag.mem")
        os.system("cp frag.mem ho_frag.mem")
        protein_size = myPersonalFunctions.file_len("ssweight")
        # print(protein_name, protein_size)
        with open("info.dat", 'w') as info:
            n = 0
            nn = 0
            with open('frag1.mem', 'w') as w:
                with open('ha_frag.mem', 'r') as f:
                    for i in range(4):
                        w.write(next(f))
                    ha_lines = f.readlines()
                # print(ha_lines)
                with open('ho_frag.mem', 'r') as f:
                    for i in range(4):
                        next(f)
                    ho_lines = f.readlines()
                for i in range(1, protein_size-7):
                    fragGroup = []
                    for line in ho_lines:
                        _, loc, _, fragLens, _ = line.split(" ")
                        locEnd = int(loc) + int(fragLens)
                        # print(loc, fragLens, locEnd)
                        if(i <= int(loc) and i+9 >= locEnd):
                            fragGroup += [line]
                    groupLens = len(fragGroup)
                    if(groupLens != 0):
                        nn += 1
                        for j in range(20):
                            w.write(fragGroup[j % groupLens])
                    else:
                        n += 1
                        for j in range(20):
                            try:
                                w.write(ha_lines[(i-1)*20+j])
                            except:
                                print(i*20+j)
                info.write("Use "+str(nn)+" HO and "+str(n)+" HA\n")
        os.chdir("../..")
if(args.jan11):
    folder_list = ["T0766"]
    # folder_list = ["T0792"]
    os.system("mkdir -p aawsemJan11")
    os.chdir("aawsemJan11")
    for protein_name in folder_list:
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        # os.system("cp  . -r".format(protein_name))
        os.system("mkdir -p simulation")
        os.chdir("simulation")
        for i in range(20):
            my_from = "../../../aawsemDec25/{0}/simulation/".format(protein_name)+str(i)
            my_to = "."
            cmd = "rsync -a --exclude='dump.lammpstrj' --exclude='slurm-*' --exclude='movie*' --exclude='q*' {} {}".format(my_from, my_to)
            print(cmd)
            os.system(cmd)
            os.chdir(str(i))
            os.system("sed -i '/read_data/c\\read_restart restart.8000000' {}.in".format(protein_name))
            os.system("sed -i 's/600/500/g' *.in")  # only apply to protein less than 200 residues, and only one .in file
            os.system("sbatch run.slurm")
            os.chdir("..")
        os.chdir("../..")

if(args.jan10):
    folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    os.system("mkdir -p aawsemJan10")
    os.chdir("aawsemJan10")
    for protein_name in folder_list:
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        # os.system("cp  . -r".format(protein_name))
        os.system("mkdir -p simulation")
        os.chdir("simulation")
        for i in range(20):
            my_from = "../../../aawsemJan08/{0}/simulation/".format(protein_name)+str(i)
            my_to = "."
            cmd = "rsync -a --exclude='dump.lammpstrj' --exclude='slurm-*' --exclude='movie*' --exclude='q*' {} {}".format(my_from, my_to)
            print(cmd)
            os.system(cmd)
            os.chdir(str(i))
            os.system("sed -i '/read_data/c\\read_restart restart.8000000' {}.in".format(protein_name))
            os.system("sed -i 's/600/500/g' *.in")  # only apply to protein less than 200 residues, and only one .in file
            os.system("sbatch run.slurm")
            os.chdir("..")
        os.chdir("../..")

if(args.jan07):
    folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    os.system("mkdir -p aawsemJan07")
    os.chdir("aawsemJan07")
    for protein_name in folder_list:
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        os.system("cp -r ../../aawsemDec25/{0}/{0} .".format(protein_name))
        os.chdir(protein_name)
        os.system("cp ~/Research/CASP11/HA_Mem/{}/frag.mem ha_frag_old.mem".format(protein_name))
        os.system("sed 's/\/work\/pw8\/mc70\/script/\/home\/wl45/g' ha_frag_old.mem > ha_frag.mem")
        os.system("cp frag.mem ho_frag.mem")
        protein_size = myPersonalFunctions.file_len("ssweight")
        # print(protein_name, protein_size)
        with open("info.dat", 'w') as info:
            n = 0
            nn = 0
            with open('frag1.mem', 'w') as w:
                with open('ha_frag.mem', 'r') as f:
                    for i in range(4):
                        w.write(next(f))
                    ha_lines = f.readlines()
                # print(ha_lines)
                with open('ho_frag.mem', 'r') as f:
                    for i in range(4):
                        next(f)
                    ho_lines = f.readlines()
                for i in range(1, protein_size-7):
                    fragGroup = []
                    for line in ho_lines:
                        _, loc, _, fragLens, _ = line.split(" ")
                        locEnd = int(loc) + int(fragLens)
                        # print(loc, fragLens, locEnd)
                        if(i <= int(loc) and i+9 >= locEnd):
                            fragGroup += [line]
                    groupLens = len(fragGroup)
                    if(groupLens != 0):
                        nn += 1
                        for j in range(20):
                            w.write(fragGroup[j % groupLens])
                    else:
                        n += 1
                        for j in range(20):
                            try:
                                w.write(ha_lines[(i-1)*20+j])
                            except:
                                print(i*20+j)
                info.write("Use "+str(nn)+" HO and "+str(n)+" HA\n")
        os.chdir("../..")

if(args.jan03):
    folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    os.system("mkdir -p aawsemJan03")
    os.chdir("aawsemJan03")
    for protein_name in folder_list:
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        # os.system("cp  . -r".format(protein_name))
        os.system("mkdir -p simulation")
        os.chdir("simulation")
        for i in range(20):
            my_from = "../../../aawsemDec25/{0}/simulation/".format(protein_name)+str(i)
            my_to = "."
            cmd = "rsync -a --exclude='dump.lammpstrj' --exclude='slurm-*' --exclude='movie*' --exclude='q*' {} {}".format(my_from, my_to)
            print(cmd)
            os.system(cmd)
            os.chdir(str(i))
            os.system("sed -i '/read_data/c\\read_restart restart.8000000' {}.in".format(protein_name))
            os.system("sed -i 's/600/500/g' *.in")  # only apply to protein less than 200 residues, and only one .in file
            os.system("sbatch run.slurm")
            os.chdir("..")
        os.chdir("../..")
if(args.dec25):
    folder_list = ["T0792", "T0815", "T0778", "T0766", "T0782", "T0833", "T0844", "T0842", "T0846", "T0803"]
    # folder_list = ["T0815", "T0778"]
    # folder_list = ["T0766"]
    name = "aawsemFeb15"
    os.system("mkdir -p {}".format(name))
    os.chdir("{}".format(name))
    for protein_name in folder_list:
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        os.system("mkdir -p "+protein_name)
        os.chdir(protein_name)
        os.system("cp ~/opt/AAWSEM/"+protein_name+".fasta .")
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
