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
import glob
from time import sleep
import fileinput
import numpy as np
from small_script.variable_test import variable_test
# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
# import re
# numbers = re.compile(r'(\d+)')
# def numericalSort(value):
#     parts = numbers.split(value)
#     parts[1::2] = map(int, parts[1::2])
#     return parts
# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="This is my playground for current project")
parser.add_argument("-r", "--run", help="test mode",
                    action="store_true")
parser.add_argument("-s", "--see", help="test mode",
                    action="store_true")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=0)
args = parser.parse_args()
if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

base_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''

if args.mode == 18:
    rg_list = [0, 0.1, 0.2, 1]
    memb_k_list = [0, 1, 2, 4]
    pressure_list = [0, 0.1, 0.2, 0.4, 0.8, 1, 2]
    # rg_list = [0.1]
    # memb_k_list = [1]
    # pressure_list = [0.1, 1]
    variable_test(rg_list=rg_list, memb_k_list=memb_k_list, pressure_list=pressure_list, repeat=2)

if args.mode == 17:
    # protocol_list = ["er", "awsemer", "frag", "raptor"]
    protocol_list = ["awsemer", "frag"]
    protein_list = ["1occ"]
    for protein in protein_list:
        for protocol in protocol_list:
            print("Work on protein: {}, protocol: {}".format(protein, protocol))
            if protocol == "raptor":
                do("cp ~/opt/gremlin/protein/1occ/raptor/go_rnativeC* {}/".format(protein))
            else:
                do("cp ~/opt/gremlin/protein/1occ/gremlin/go_rnativeC* {}/".format(protein))
            do("mkdir -p {}".format(protocol))
            do("cp -r {} {}/".format(protein, protocol))
            cd(protocol)
            cd(protein)
            fileName = "{}_multi.in".format(protein)
            if protocol == "raptor":
                backbone_file = "fix_backbone_coeff_er.data"
                do("cp ~/opt/gremlin/protein/{}/raptor/go_rnativeC* .".format(protein))
            else:
                backbone_file = "fix_backbone_coeff_{}.data".format(protocol)
                do("cp ~/opt/gremlin/protein/{}/gremlin/go_rnativeC* .".format(protein))
            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    tmp = line
                    tmp = tmp.replace("fix_backbone_coeff_er.data", backbone_file)
                    print(tmp, end='')
            cd("..")
            do("run.py -m 0 -n 20 {}".format(protein))
            cd("..")

if args.mode == 16:
    rg_list = [0, 0.1, 0.2, 0.4, 0.5, 1, 2, 4]
    variable_test(rg_list=rg_list, repeat=1, commons=True)

if(args.mode == 15):
    print("create directory_list")
    with open("directory_list", "w") as f:
        for i in range(40):
            # print(os.getcwd())
            location = os.getcwd() + "/../"
            f.write(location+str(i)+"/0\n")
    do("cp ../../2xov/2xov.pdb .")
    do("python2 ~/opt/small_script/CalcLocalDistanceStats.py 2xov directory_list out")
if(args.mode == 14):
    print("Extract qw and distance info.")
    for i in range(100):
        cd(str(i))
        cd("0")
        do("awk '{print $2}' wham.dat |  sed 's/,$//' | sed 1d > qw.dat")
        do("awk '{print $2}' addforce.dat |  sed 's/,$//' | sed 1d > distance.dat")
        cd("../..")

if args.mode == 13:
    rg_list = [0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
    memb_k_list = [0, 1, 2, 4, 8]
    variable_test(rg_list=rg_list, memb_k_list=memb_k_list)
if args.mode == 12:
    rg_list = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
    variable_test(rg_list=rg_list)
if args.mode == 11:
    zim_type_list = ["aug04", "aug26"]
    membrane_width_list = [30, 28.8]
    for zim in zim_type_list:
        for width in membrane_width_list:
            folder = "zim_{}_width_{}".format(zim, width)
            do("mkdir -p {}".format(folder))
            cd(folder)
            do("cp -r ../2xov .")
            cd("2xov")
            fixFile = "fix_backbone_coeff_single.data"
            with fileinput.FileInput(fixFile, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace("WIDTH", str(width)), end='')
            do("cp zim_{} zim".format(zim))
            cd("..")
            do("run.py -n 2 2xov")
            cd("..")
if args.mode == 10:
    distance_list = np.linspace(166, 180, 15)
    for distance in distance_list:
        folder = "dis_{}".format(distance)
        cd(folder)
        do("sbatch run_0.slurm")
        cd("..")
if args.mode == 9:
    cmd = "python3 ~/opt/small_script/find_distance.py"
    run_slurm = base_slurm.format(cmd)
    folder_list = ['force_0.045']
    print(folder_list)
    for folder in folder_list:
        cd(folder)
        cd("simulation")
        run_list = glob.glob("*")
        for run in run_list:
            cd(run)
            cd("0")
            with open("find_distance.slurm", "w") as r:
                r.write(run_slurm)
            do("sbatch find_distance.slurm")
            cd("../..")
        cd("../..")

if args.mode == 8:
    cmd = "gg.py -m 8"
    run_slurm = base_slurm.format(cmd)

    # folder_list = glob.glob("force_*")
    # folder_list = ['force_0.08', 'force_0.03', 'force_0.0']
    folder_list = ['force_0.055']
    # folder_list = ['force_0.07', 'force_0.02', 'force_0.045']
    # folder_list = ['force_0.06', 'force_0.04']
    print(folder_list)
    for folder in folder_list:
        cd(folder)
        cd("simulation")
        run_list = glob.glob("*")
        for run in run_list:
            cd(run)
            cd("0")
            with open("compute_angle.slurm", "w") as r:
                r.write(run_slurm)
            do("sbatch compute_angle.slurm")
            cd("../..")
        cd("../..")

if args.mode == 7:
    for i in range(80):
        do("mv {0} ../../../new_force_ramp/memb_0_force_ramp_rg_0_new/simulation/{1}".format(i,i+90))
if args.mode == 6:
    force_list = [0.55, 0.6, 0.65]
    # force_list = [0.25, 0.35, 0.4, 0.45]
    # force_list = [0.15, 0.2]
    for force in force_list:
        do("mkdir force_{}".format(force))
        do("cp -r 2xov force_{}/".format(force))
        cd("force_{}".format(force))
        with fileinput.FileInput("2xov/2xov_multi.in", inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("MY_FORCE", str(force)), end='')
        do("run.py -n 10 2xov/")
        cd("..")


if args.mode == 5:
    # cd("start_misfolded")
    distance_list = np.linspace(0, 30, 16)
    for dis in distance_list:
        do("mkdir -p dis_{}".format(dis))
        do("cp -r ../2xov/ dis_{}".format(dis))
        do("cp ../../freeEnergy/go_model_start_unfolded/simulation/dis_{0}/restart.25000000 dis_{0}/2xov/".format(dis))
        cd("dis_{}".format(dis))
        do("run.py -n 10 2xov/")
        cd("..")

if args.mode == 4:
    do("rm data")
    for i in range(100):
        do("cat dis_{}/0/data >> data.dat".format(i))
    do("awk '{print $1}' data.dat  > e.dat")
    do("awk '{print $2}' data.dat  > p.dat")
    do("awk '{print $3}' data.dat  > qw.dat")

if args.mode == 1:
    cd("simulation")
    do("pulling_prepare.py")
    cd("..")
    do("mkdir freeEnergy")
    cd("freeEnergy")
    do("make_metadata.py -k 0.05 -t 300")
    do("pulling_analysis.py -m 3 -p 2")

if args.mode == 2:
    print("80 bins.")
    # cd("simulation")
    # do("pulling_prepare.py")
    # cd("..")
    do("mkdir more_bin")
    cd("more_bin")
    do("make_metadata.py -k 0.05 -t 600")
    do("pulling_analysis.py -m 3 -p 1")

if args.mode == 3:
    # cd("simulation")
    # do("pulling_prepare.py")
    # cd("..")
    do("mkdir -p only_less_than_100")
    cd("only_less_than_100")
    do("make_metadata.py -k 0.05 -t 300 -m 2")
    do("pulling_analysis.py -m 3 -p 1")
    # for i in range(52, 70):
    #     do("mv {}/{} .".format(i, i-40))
    #     do("mv {} {}".format(i, i-20))
    # for i in range(50):
    #     do("mv force_0.8_2/{} force_0.8/{}".format(i, i+50))
        # do("mv half_contact_force_0.8_memb1_rg1_2/{} half_contact_force_0.8_memb1_rg1/{}".format(i, i+20))


if(args.run):
    print("Hello World")

    name = "T0833"
    n = 21
    do("mkdir "+name)
    cd(name)
    for i in range(1, n):
        do("mkdir -p job.{}".format(i))
        do("cp ../preparation_files/myjob_nots.slurm job.{}".format(i))
        do("cp ../preparation_files/loopsubmit.bash job.{}".format(i))
        do("cp -r ../preparation_files/{0}_runpackage job.{1}/runpackage".format(name, i))
        do("cp ../preparation_files/{1}_tpr/run.{0}.tpr job.{0}/runpackage/run.tpr".format(i, name))
    for i in range(1, n):
        cd("job.{}".format(i))
        fileName = "myjob_nots.slurm"
        with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("T0766", name), end='')

        do("bash loopsubmit.bash")
        cd("..")

    # for i in range(1, 6):
    #     do("mkdir job.{}".format(i))
    #     do("cp myjob_nots.slurm job.{}".format(i))
    #     do("cp loopsubmit.bash job.{}".format(i))
    #     do("cp -r runpackage job.{}".format(i))
    #     do("cp run.{0}.tpr job.{0}/runpackage/run.tpr".format(i))
    # for i in range(1, 6):
    #     cd("job.{}".format(i))
    #     fileName = "myjob_nots.slurm"
    #     name = "T0833"
    #     with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
    #         for line in file:
    #             print(line.replace("T0766", name), end='')
    #     do("bash loopsubmit.bash")
    #     cd("..")

    # force_list = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
    # for force in force_list:
    #     folder = "1d_force_" + str(force)
    #     do("mkdir -p " + folder)
    #     cd(folder)
    #     do("cp ../metadatafile .")
    #     do("~/bin/python3/bin/python3 ~/opt/pulling_analysis.py -f -m 5 --force {}".format(force))
    #     do("sbatch freeEnergy.slurm")
    #     cd("..")
    # do("rm -r test")
    # do("cp -r 2xov test")
    # cd("test")
    # do("test_run.py test.in")

if(args.see):
    do("head test/0/addforce.dat")
def fix_error_run():
    n = args.number
    for i in range(n):
        os.chdir(str(i))
        os.system("cp ~/opt/pulling/qnqc.slurm .")
        os.system("sbatch qnqc.slurm")
        os.chdir("..")
    # os.system("grep 'srun: error: Application launch failed: Socket timed out on send/recv operation' . -r | cut -d':' -f1 | rev | cut -d"/" -f2- | rev > list")
    # array = []
    # cwd = os.getcwd()
    # print(cwd)
    # with open('list', 'r') as ins:
    #     for line in ins:
    #         target = line.strip('\n')
    #         array.append(target)
    # for i in array:
    #     os.chdir(i)
    #     os.system("rm slurm-*")
    #     os.system("sbatch rerun.slurm")
    #     sleep(0.5)  # Time in seconds.
    #     os.system("sbatch qnqc.slurm")
    #     sleep(0.5)  # Time in seconds.
    #     # os.system("pwd")
    #     os.chdir(cwd)

    # os.system("cut -d'/' -f2 list >")
# if(args.fix):
#     fix_error_run()


def rerun():
    n = args.number
    for i in range(n):
        os.system("cp -r {0} rerun_{0}".format(str(i)))
        # source = "~/opt/gagb/gagb_constant200_rerun.in"
        # target = " 2lhc.in"
        source = "~/opt/pulling/2xov_force_load_dis.in"
        target = " 2xov.in"
        os.system("cp "+source + target)
        os.system("cp "+target+" rerun_{}/".format(str(i)))
        os.chdir("rerun_"+str(i))
        os.system("rm slurm*")
        os.system("sbatch run.slurm")
        os.chdir("..")
# if(args.rerun):
#     rerun()


def continue_run():
    n = 10
    os.chdir("simulation")
    for i in range(n):
        os.system("cp -r {0} continue_{0}".format(str(i)))
        # source = "~/opt/gagb/gagb_constant200_rerun.in"
        # target = " 2lhc.in"
        source = "~/opt/pulling/2xov_continue_run.in"
        target = " 2xov.in"
        os.system("cp "+source + target)
        os.system("cp "+target+" continue_{}/".format(str(i)))
        os.chdir("continue_"+str(i))
        os.system("rm slurm*")
        os.system(  # replace RANDOM with a radnom number
            "sed -i.bak 's/RANDOM/'" +
            str(randint(1, 10**6)) +
            "'/g' 2xov.in")
        os.system("sbatch run.slurm")
        os.chdir("..")
# if(args.go):
#     continue_run()
# parser = argparse.ArgumentParser(
#         description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
# # parser.add_argument("template", help="the name of template file")
# args = parser.parse_args()
# # protein_name = args.template.split('_', 1)[-1].strip('/')
# protein_name = args.protein.strip('/')
    # name = "ga_2m"
## -------------Pulling--------
# os.system("cp ~/opt/small_script/springForce.plt .")
# os.system("cp ~/opt/small_script/springForce_smooth.plt .")
# os.system("gnuplot springForce.plt")
# os.system("gnuplot springForce_smooth.plt")
# os.system("open springForce.pdf")
# os.system("open springForce_smooth.pdf")
# SpringConstant_list = [3e-05, 5e-05, 1e-06, 3e-06, 5e-06, 1e-05, 1e-07]
# # SpringConstant_list = [3e-06, 5e-06]
# for SpringConstant in SpringConstant_list:
#     name = "spring"+str(SpringConstant)
#     os.system("mkdir "+name)
#     os.chdir(name)
#     os.system("cp -r ../2xov/ .")
#     os.system("cp ../variables.dat .")
#     os.chdir("2xov")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#         "sed -i.bak 's/SpringForce/'" +
#         str(SpringConstant) +
#         "'/g' "+protein_name+".in")
#     os.chdir("..")
#     os.system("run.py 2xov/ -s 8 -n 2")
#     os.chdir("..")

# number_of_run_list = [2, 4, 8, 16]
# for n in number_of_run_list:
#     name = "ga_"+str(n)+"m"
#     # os.system("mkdir "+name)
#     os.system("cp -r 2lhd.pdb "+name)
#
#     # os.system("cp -r 2lhc variables.dat "+name)
#     os.chdir(name)
#     for i in range(20):
#         os.chdir("analysis/"+str(i))
#         os.system("cp ../../2lhd.pdb .")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhd.pdb dump.lammpstrj q_gb.dat")
#         os.system("python2 ~/opt/script/CalcQValue.py 2lhc.pdb dump.lammpstrj q_ga.dat")
#         os.system("cp ~/opt/small_script/qw_gagb.plt .")
#         os.system("gnuplot qw_gagb.plt")
#         os.system("mv qw_gagb.pdf ../../results/qw_gagb_{0}.pdf".format(str(i)))
#         os.chdir("../..")
#     os.chdir("..")
#
# for n in number_of_run_list:
#     name = "ga_"+str(n)+"m"
#     # os.system("mkdir "+name)
#     os.system("cp -r 2lhd.pdb "+name)
#
#     # os.system("cp -r 2lhc variables.dat "+name)
#     os.chdir(name)
#     for i in range(20):
#         os.chdir("analysis/"+str(i))
#         os.system("paste q_ga.dat q_gb.dat > q_gagb.dat")
#         os.system("cp ~/opt/small_script/qw_ga-gb.plt .")
#         os.system("gnuplot qw_ga-gb.plt")
#         os.system("mv qw_ga-gb.pdf ../../results/qw_ga-gb_{0}.pdf".format(str(i)))
#         os.chdir("../..")
#     os.system("cp ~/opt/small_script/qw_ga_all.plt .")
#     os.system("gnuplot qw_ga_all.plt")
#     os.system("cp ~/opt/small_script/qw_gb_all.plt .")
#     os.system("gnuplot qw_gb_all.plt")
#     os.system("cp ~/opt/small_script/qw_diff_all.plt .")
#     os.system("gnuplot qw_diff_all.plt")
#     os.chdir("..")

# simulation_steps = 4 * 10**6
# warm_up_steps = 10 * 10**5
#
# seed(datetime.now())
# n= 20
# vmd = "/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command"
#
# os.system("BuildAllAtomsFromLammps.py dump.lammpstrj movie")
# os.system("cp ~/opt/plot_scripts/2xov_movie.tcl .")
# os.system(vmd+" -e 2xov_movie.tcl ")
# os.system("mkdir -p MyResults")
# for i in range(n):
#     print(i)
#     os.chdir("analysis/"+str(i))
#     os.system("cp ~/opt/plot_scripts/2xov_movie_screenshot.tcl .")
#     os.system(vmd+" -e 2xov_movie_screenshot.tcl")
#     os.system("cp frame1000.tga ../../MyResults/frame"+str(i)+"_1000.tga")
#     #os.system("cp frame450.tga ../Results/frame"+folder_name+"_450.tga")
#     # os.system("movie.py "+protein_name)
#     os.chdir("../..")
#     # analysis

# folder_name = ""
# result_folder = "WeiLu_Aug_07"

# protein_list = ['T089', 'T120', 'T251', 'TOP7', '1UBQ']
# sublist = ['']
# # sublist = ['_ha', '_he']
# # sublist = ['_lp', '_he_lp']
# # folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_name = protein+sub
#         os.chdir(folder_name)
#         os.chdir("best_2nd")
#         os.system("pymol ~/opt/plot_scripts/align.pml > matrix.dat")
#         os.system("head -n 70 matrix.dat | tail -n 20 > cealign_matrix.dat")
#         # for i in range(19, -1, -1):
#         #     os.system("mv {}.pdb {}.pdb".format(i, i+1))
#         os.chdir("../..")
    # os.chdir(protein)
    # os.chdir("best_1st")
    # os.system("python3 ~/opt/small_script/cross_q.py")
    # os.chdir("..")
    # os.chdir("best_2nd")
    # os.system("python3 ~/opt/small_script/cross_q.py")
    # os.chdir("..")
    # os.chdir("..")
# n = 3
# for i in range(n):
#     # simulation set up
#     folder_name = str(i)
#     os.system("mkdir -p "+folder_name)
#     os.system("cp -r "+args.protein+"* "+folder_name)
#     os.chdir(folder_name)
#     os.system("cp ../../helix_less/simulation/"+str(i)+"/restart.4000000 .")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#         "sed -i.bak 's/WARM_UP_STEPS/'" +
#         str(warm_up_steps) +
#         "'/g' "+protein_name+".in")
#     os.system(  # replace RANDOM with a radnom number
#             "sed -i.bak 's/RANDOM/'" +
#             str(randint(1, 10**6)) +
#             "'/g' "+protein_name+".in")
#     os.system(  # replace SIMULATION_STEPS with specific steps
#             "sed -i.bak 's/SIMULATION_STEPS/'" +
#             str(simulation_steps) +
#             "'/g' "+protein_name+".in")
# # if(platform.system() == 'Darwin'):
# #     os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
# #     < "+protein_name+".in")
#     if(platform.system() == 'Darwin'):
#         os.system("/Users/weilu/Documents/lammps-9Oct12_modified/src/lmp_serial \
#         < "+protein_name+".in")
#     elif(platform.system() == 'Linux'):
#         os.system("cp ~/opt/run.slurm .")
#         os.system(  # replace PROTEIN with pdb name
#                 "sed -i.bak 's/PROTEIN/'" +
#                 protein_name +
#                 "'/g' run.slurm")
#         os.system("sbatch run.slurm")
#     else:
#         print("system unkown")
#     os.chdir("..")
# exit(1)

# w_helix_list = [0.1, 0.5, 1, 1.5]
# m_helix_list = [0.1, 0.5, 1, 1.5]
#
# for i in range(len(w_helix_list)):
#     w = w_helix_list[i]
#     for j in range(len(m_helix_list)):
#
#         # m = m_helix_list[j]
#         folder_name = str(i)+"_"+str(j)
#         # os.system("cd "folder_name)
#         os.chdir(folder_name)
#         # os.system("analysis.py 2xov/")
#         # os.system("echo "+folder_name+" >> ../all")
#         os.system("sort -k 3 analysis/list_of_max_q > ../data/"+folder_name)
#         os.chdir("..")
#         # os.system("mkdir "+folder_name)
#         # os.chdir(folder_name)
#         # os.system("cp -r ../2xov .")
#         # os.chdir("2xov")
#         # os.system(
#         #         "sed -i.bak 's/W_HELIX/'" +
#         #         str(w) +
#         #         "'/g' fix_backbone_coeff.data")
#         # os.system(
#         #         "sed -i.bak 's/M_HELIX/'" +
#         #         str(m) +
#         #         "'/g' fix_backbone_coeff.data")
#         # os.chdir("..")
#         # os.system("run.py 2xov/ -n 5")

# os.system("cp ~/opt/gg.py this_gg.py")
# for i in range(5):
#     os.system("mkdir "+str(i))
#     os.chdir(str(i))
#     os.system("cp -r ../2xov/ .")
#     os.system("cp ../../2xov_strong_single_memory_600to500/simulation/"+str(i)+"/restart.2000000 2xov/")
#     os.system("run.py -s 4 -n 2 2xov/")
#     os.chdir("..")

# # rama_list = [6, 8, 16]
# # rama_list = [4]
# melt_t_list = [400, 500, 600]
# for variable in melt_t_list:
#     folder_name = str(variable)
#     os.system("mkdir "+folder_name)
#     os.chdir(folder_name)
#     os.system("cp -r ../1qjp .")
#     os.chdir("1qjp")
#     os.system(
#             "sed -i.bak 's/MELTT/'" +
#             str(variable) +
#             "'/g' 1qjp.in")
#     os.chdir("..")
#     # os.system("pwd")
#     os.system("run.py 1qjp/ -n 5 -s 5")
#     os.chdir("..")
# os.system("cp ~/opt/gg.py this_gg.py")
#
# exec (open("config.py").read())
# n = number_of_run
# steps = simulation_steps
#
# protein_name = args.protein.strip('/')
#
# temp = 400
# folder_name = "{}_t{}_q100_test11".format(protein_name, str(temp))
# print("all going to "+folder_name)
# os.system("mkdir -p "+folder_name)
# os.system("rm -f "+folder_name + "/*")
# command = 'cat simulation/{}/%d/wham11 \
# >> {}/all_wham.dat'.format(temp, folder_name)
# # cal rmsd
# os.chdir("simulation/"+str(temp))
# for i in range(n):
#     os.chdir(str(i))
#     os.system("awk '{print>\"file1\"(NR>(n/2)?2:1)}' n=\"$(wc -l <file1)\" file1")
#     os.system("cat file11 >> ../../../"+folder_name+"/rmsd_total")
#     # os.system("sed 1d wham.dat > wham1d.dat")
#     os.system("awk '{print>\"wham1\"(NR>(n/2)?2:1)}' n=\"$(wc -l <wham1)\" wham1")
#     os.chdir("..")
# os.chdir("../..")
# for i in range(n):
#     cmd = command % i
#     os.system(cmd)
# os.chdir(folder_name)
# os.system("awk '{print $2}' all_wham.dat > Qw_total")
# os.system("awk '{print $3}' all_wham.dat > rg_total")
# os.system("awk '{print $4}' all_wham.dat > p_total")
# os.system("awk '{print $5}' all_wham.dat > tc_total")
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# os.system("cp ~/opt/wham_analysis/*.m .")
# os.chdir("..")
# os.system("~/opt/script/wham/fused_calc_cv.sc {} top7 50 400 350 450 5 50 100 0 0.98".format(folder_name))
#
#
# folder_name = "{}_t{}_q100_test12".format(protein_name, str(temp))
# print("all going to "+folder_name)
# os.system("mkdir -p "+folder_name)
# os.system("rm -f "+folder_name + "/*")
# command = 'cat simulation/{}/%d/wham12 \
# >> {}/all_wham.dat'.format(temp, folder_name)
# # cal rmsd
# os.chdir("simulation/"+str(temp))
# for i in range(n):
#     os.chdir(str(i))
#     os.system("cat file12 >> ../../../"+folder_name+"/rmsd_total")
#     os.chdir("..")
# os.chdir("../..")
# for i in range(n):
#     cmd = command % i
#     os.system(cmd)
# os.chdir(folder_name)
# os.system("awk '{print $2}' all_wham.dat > Qw_total")
# os.system("awk '{print $3}' all_wham.dat > rg_total")
# os.system("awk '{print $4}' all_wham.dat > p_total")
# os.system("awk '{print $5}' all_wham.dat > tc_total")
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# os.system("cp ~/opt/wham_analysis/*.m .")
# os.chdir("..")
#
#
#
# os.system("~/opt/script/wham/fused_calc_cv.sc {} top7 50 400 350 450 5 50 100 0 0.98".format(folder_name))

#
# result_folder = "WeiLu_Aug_07"
# os.system("mkdir -p "+result_folder)
# protein_list = ['T089', 'T120', 'T251', 'top7', '1UBQ']
# # sublist = ['_ha', '_he']
# sublist = ['_lp', '_he_lp']
# folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_list += [protein+sub]
# print(folder_list)
# # exit(1)
# # awk '{print>'file'(NR>(n/2)?2:1)}' n='$(wc -l <test)' test
# for folder in folder_list:
#     print(folder)
#     os.chdir(folder)
#     exec (open("config.py").read())
#     n = number_of_run
#     steps = simulation_steps
#     os.system("mkdir -p ../{}/".format(result_folder)+folder+"/best_q")
#     os.system("sort analysis/list_of_max_q > ../{}/q_".format(result_folder)+folder+".dat")
#     for i in range(n):
#         # move
#         os.chdir("analysis/"+str(i))
#         os.system("cp chosen.pdb ../../../{}/".format(result_folder) + folder+"/best_q/"+str(i)+".pdb")
#         os.chdir("../..")
#     os.chdir("..")


# result_folder = "WeiLu_Aug_07"
# os.system("mkdir -p "+result_folder)
# protein_list = ['T089', 'T120', 'T251', 'top7', '1UBQ']
# # sublist = ['_ha', '_he']
# sublist = ['_lp', '_he_lp']
# folder_list = []
# for protein in protein_list:
#     for sub in sublist:
#         folder_list += [protein+sub]
# print(folder_list)
# # exit(1)
#
# for folder in folder_list:
#     print(folder)
#     os.chdir(folder)
#     exec (open("config.py").read())
#     n = number_of_run
#     steps = simulation_steps
#     os.system("mkdir -p ../{}/".format(result_folder)+folder+"/best_q")
#     os.system("sort analysis/list_of_max_q > ../{}/q_".format(result_folder)+folder+".dat")
#     for i in range(n):
#         # move
#         os.chdir("analysis/"+str(i))
#         os.system("cp chosen.pdb ../../../{}/".format(result_folder) + folder+"/best_q/"+str(i)+".pdb")
#         os.chdir("../..")
#     os.chdir("..")
