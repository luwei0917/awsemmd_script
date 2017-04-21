#!/usr/bin/env python2
import os, sys, re

# Scripts and programs
python_exec = 'python2'
build_pdb_script = '~/opt/script/BuildAllAtomsFromLammps.py'

# Auxillary functions
def build_pdb(dump_file, snapshot_index, structure_index):
    path, file_name = os.path.split(dump_file)
    print path
    print "%s %s %s structure%s %s -seq %s/../2xov.seq" % (python_exec, build_pdb_script, dump_file, structure_index, snapshot_index, path)
    os.system("%s %s %s structure%s %s -seq ~/opt/pulling/2xov.seq" % (python_exec, build_pdb_script, dump_file, structure_index, snapshot_index))
# render Tachyon select.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o select.tga -res 2000 2000
# Lists
data_array = []
files_array = []
conditions = []
condition_signs = []

# File names and parameters
dump_file_name = "dump.lammpstrj"
output_file_name = "pick_structures_all.dat"
structure_output_file_name = "pick_structures.dat"
structure_index = 1
structure_stride = 5050
max_pdbs_to_build = 20
pdb_index = 0
found_pdb_index = 0

if len(sys.argv) < 4:
    print "Usage: python pick_structures.py metadata output_directory cond1 (cond2 cond3 ...)"
    sys.exit()

# Read command line arguments
metadata_file = sys.argv[1]
output_directory = os.path.abspath(sys.argv[2])

if len(sys.argv) > 3:
    for j in range(3,len(sys.argv)):
        if "gt" in sys.argv[j]:
            condition_signs.append("+")
            condition = sys.argv[j].split("gt")
        if "lt" in sys.argv[j]:
            condition_signs.append("-")
            condition = sys.argv[j].split("lt")
        conditions.append(condition)

# print(conditions)
# Load all data into array
for line in open(metadata_file, "r"):
    line = line.split()
    file_name = line[0]
    if file_name == "#": continue
    files_array.append(file_name)
    data_array.append([])
    line_index = 0
    for data in open(file_name, "r"):
        data = data.strip("\n").split(",")
        # print(data)
        if data[0] == "#timestep": continue
        if data[0] == "Steps": continue
        data_array[files_array.index(file_name)].append([])
        for i in range(len(data)):
            data_array[files_array.index(file_name)][line_index].append(float(data[i]))
        line_index += 1
    # print(len(data_array))
# loop over data and output those points that satisfy all conditions
if not os.path.exists(output_directory):
    os.mkdir(output_directory)
os.chdir(output_directory)
output_file = open(output_file_name, "w")
structure_output_file = open(structure_output_file_name, "w")
vmd_out = open("vmd.tcl", "w")

for data_file in files_array:
    new_traj = True
    path_name, data_file_name = os.path.split(data_file)
    file_index = files_array.index(data_file)
    data_index = 0
    for data_point in data_array[file_index]:
        bad_condition = False
        # If all conditions are satisfied, print out the data
        for i in range(len(conditions)):
            if condition_signs[i] == "+":
                if not data_point[int(conditions[i][0])-1] > float(conditions[i][1]): bad_condition = True
            elif condition_signs[i] == "-":
                if not data_point[int(conditions[i][0])-1] < float(conditions[i][1]): bad_condition = True
            else:
                print "Bad condition argument."
                sys.exit()
        if not bad_condition:
            # print int(data_array[file_index].index(data_point)+1), found_pdb_index
            found_pdb_index += 1
            dump_file = data_file.replace(data_file_name,dump_file_name)
            output_file.write("%d %s %s\n" % (structure_index, dump_file, str(data_point).replace(',','').replace('[','').replace(']','')))
            # if pdb_index < max_pdbs_to_build and (found_pdb_index % structure_stride == 0 or new_traj is True):
            if pdb_index < max_pdbs_to_build and found_pdb_index % structure_stride == 0:
                structure_output_file.write("%d %s %s\n" % (structure_index, dump_file, str(data_point).replace(',','').replace('[','').replace(']','')))
                build_pdb(dump_file, int(data_array[file_index].index(data_point)+1), structure_index)
                vmd_out.write("mol new structure%s.pdb\n" % structure_index)
                vmd_out.write("mol modcolor 0 [molinfo top] Index\n")
                vmd_out.write("mol modstyle 0 [molinfo top] NewCartoon 0.300000 10.000000 4.100000 0\n")
                pdb_index += 1
                new_traj = False
            structure_index += 1

        data_index += 1

output_file.close()
vmd_out.write("color Display Background white\n")
vmd_out.close()
