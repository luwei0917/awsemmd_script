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
from time import sleep
import fileinput

if(platform.system() == 'Darwin'):  # Mac system (local machine)
    OPENAWSEM_LOCATION = "/Users/weilu/openmmawsem/"
elif(platform.system() == 'Linux'):
    OPENAWSEM_LOCATION = '/projects/pw8/wl45/openmmawsem/'
else:
    print("system unknown")
sys.path.insert(0, OPENAWSEM_LOCATION)
from openmmawsem import *


# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulations")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("--name", default="simulation", help="Name of the simulation")
parser.add_argument("-n", "--number", type=int, default=1,
                    help="# of simulation run, default: 1")
parser.add_argument("--restart", type=int, default=0,
                    help="start from? default: 0")
parser.add_argument("--runs", type=int, default=1,
                    help="then do how many runs?, default: 1")
parser.add_argument("-s", "--steps", type=int, default=50,
                    help="How many steps in unit of hundred thousand(not millions),\
                    per run, default: 50")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=4)
parser.add_argument("-i", "--inplace", action="store_true", default=False)
parser.add_argument("-f", "--force", type=float, default=1.0)
parser.add_argument("--start", default="native")
parser.add_argument("--commons", type=int, default=0)
parser.add_argument("-c", "--chain", type=str, default="A")
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

with open('commandline_args.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


if args.mode == 4:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --constraint=skylake
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/sep03/src/lmp_serial -in {}_{}.in
    '''

# SBATCH --constraint=skylake
if args.mode == 3:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --constraint=skylake
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/lammps-16Mar18/src/lmp_serial -in {}_{}.in
    '''

if args.mode == 2:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_serial -in {}_{}.in
    '''
if args.mode == 1:
    run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --constraint=opath
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_serial -in {}_{}.in
    '''

if args.commons == 1:
    run_slurm = run_slurm.replace("ctbp-common", "commons")
# run with scavenge
if args.commons == 2:
    run_slurm = run_slurm.replace("--partition=ctbp-common", "--partition=scavenge")
    run_slurm = run_slurm.replace("--account=ctbp-common", "--account=commons")
    run_slurm = run_slurm.replace("#SBATCH --constraint=skylake", "")
    run_slurm = run_slurm.replace("--time=1-00:00:00", "--time=04:00:00")


proteinName = args.protein.strip("/.")
def set_up():
    seed(datetime.now())
    with open("my_loop_submit.bash", "w") as f:
        steps = args.steps*1e5
        runs = args.runs
        restart = args.restart
        for ii in range(runs):
            i = ii + restart
            with open("run_{}.slurm".format(i), "w") as r:
                r.write(run_slurm.format(proteinName, i))
            if(i != restart):
                dependency = "--dependency=afterany:$jobid"
            else:
                dependency = ""
            f.write("jobid=`sbatch "+dependency+" run_{}.slurm".format(i) + " | tail -n 1 | awk '{print $4}'`\necho $jobid\n")
            if i == 0:
                start_from = "read_data data.{}".format(proteinName)
                if args.start == "extended":
                    start_from = "read_restart restart.extended"
                if args.start == "crystal":
                    start_from = "read_data data.crystal"
            else:
                start_from = "read_restart restart." + str(int(steps*i))
            do("cp {0}_multi.in {0}_{1}.in".format(proteinName, i))
            fileName = "{0}_{1}.in".format(proteinName, i)
            # backbone_file = "fix_backbone_coeff_{}.data".format(args.name)
            with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
                for line in file:
                    tmp = line.replace("read_data data.crystal", "START_FROM")  # remove in future.
                    tmp = tmp.replace("langevin 800 800", "langevin 300 300")  # change temp, remove in future
                    tmp = tmp.replace("langevin 800 300", "langevin 300 300")  # change temp, remove in future
                    tmp = tmp.replace("START_FROM", start_from)
                    tmp = tmp.replace("MY_FORCE", str(args.force))
                    # tmp = tmp.replace("fix_backbone_coeff_er.data", backbone_file)
                    if i == 0:
                        tmp = tmp.replace("RESET_TIME_OR_NOT", "reset_timestep	0")
                    print(tmp, end='')
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/SIMULATION_STEPS/'" +
                str(int(steps)) +
                "'/g' {}_{}.in".format(proteinName, i))
            os.system(  # replace SIMULATION_STEPS with specific steps
                "sed -i.bak 's/NUMBER/'" +
                str(int(i)) +
                "'/g' {}_{}.in".format(proteinName, i))
            os.system("mkdir -p {}".format(i))
            os.system(  # replace RANDOM with a radnom number
                "sed -i.bak 's/RANDOM/'" +
                str(randint(1, 10**6)) +
                "'/g' *.in")



def batch_run():
    if(platform.system() == 'Darwin'):
        os.system("/Users/weilu/bin/lmp_serial < "+protein_name+".in")
        # os.system("/Users/weilu/Research/Build/lammps-9Oct12_modified/src/lmp_serial \
        # < "+protein_name+".in")
    elif(platform.system() == 'Linux'):
        os.system("bash my_loop_submit.bash")
        sleep(0.2)  # Time in seconds.
    else:
        print("system unkown")

# if(args.inplace):
#     print("inplace")
#     set_up()
#     # batch_run()
#     do("~/build/brian/z_dependence/lmp_serial -in {}_multi.in".format(proteinName))
# else:
#     n = args.number
#     cwd = os.getcwd()
#     for i in range(n):
#         if args.restart == 0:
#             do("mkdir -p " + args.name)
#             do("cp -r {} {}/{}".format(proteinName, args.name, i))
#         cd(args.name + "/"+str(i))
#         set_up()
#         batch_run()
#         cd(cwd)

# print("hello world")


simulation_platform = "CPU"  # OpenCL, CUDA, CPU, or Reference
platform = Platform.getPlatformByName(simulation_platform)
platform.setPropertyDefaultValue("Threads", "1")
print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")
proteinName = pdb_id = args.protein
chain=args.chain.upper()
pdb = f"{pdb_id}.pdb"

input_pdb_filename, cleaned_pdb_filename = prepare_pdb(pdb, chain)
ensure_atom_order(input_pdb_filename)
getSeqFromCleanPdb(input_pdb_filename, chains='A')

def add_chain_to_pymol_pdb(location):
    # location = "/Users/weilu/Research/server/nov_2018/openMM/random_start/1r69.pdb"
    with open("tmp", "w") as out:
        with open(location, "r") as f:
            for line in f:
                info = list(line)
                if len(info) > 21:
                    info[21] = "A"
                out.write("".join(info))
    os.system(f"mv tmp {location}")


reporter_frequency = 4000
oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, xml_filename=OPENAWSEM_LOCATION+"awsem.xml") # k_awsem is an overall scaling factor that will affect the relevant temperature scales

# apply forces
forces = [
    oa.con_term(),
    oa.chain_term(),
    oa.chi_term(),
    oa.excl_term(),
    oa.rama_term(),
    oa.rama_proline_term(),
    oa.rama_ssweight_term(),
    oa.contact_term(),
    # oa.direct_term(),
    # oa.burial_term(),
    # oa.mediated_term(),
    oa.fragment_memory_term(frag_location_pre="./")
]
oa.addForces(forces)

# start simulation
collision_rate = 5.0 / picoseconds

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
simulation.context.setPositions(oa.pdb.positions) # set the initial positions of the atoms
# simulation.context.setVelocitiesToTemperature(300*kelvin) # set the initial velocities of the atoms according to the desired starting temperature
simulation.minimizeEnergy() # first, minimize the energy to a local minimum to reduce any large forces that might be present
simulation.reporters.append(StateDataReporter(stdout, reporter_frequency, step=True, potentialEnergy=True, temperature=True)) # output energy and temperature during simulation
simulation.reporters.append(PDBReporter("movie.pdb", reporter_frequency)) # output PDBs of simulated structures
simulation.step(int(1e5))
# simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_reporter_frequency)) # save progress during the simulation

