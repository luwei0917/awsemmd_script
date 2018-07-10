#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import myPersonalFunctions
import fileinput

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--frag", action="store_true", default=False)
parser.add_argument("--crystal", action="store_true", default=False)
parser.add_argument("--membrane", action="store_true", default=False)
parser.add_argument("--globular", action="store_true", default=False)

args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

proteinName = args.protein



# get fasta, pdb, seq file redy

do("~/opt/script/pdb2fasta.sh crystal_structure.pdb > {0}.fasta".format(proteinName))
size = myPersonalFunctions.length_from_fasta("{0}.fasta".format(proteinName))
if False:  # used for go model
    do("~/opt/fasta2pdb.py "+proteinName)
    do("python2 ~/opt/script/GetCACADistancesFile.py crystal_structure native.dat")
    do("python2 ~/opt/script/GetCACoordinatesFromPDB.py crystal_structure nativecoords.dat")
    do("cp native.dat rnative.dat")  # q bias need rnative
    # do("python2 ~/opt/Pdb2Gro.py {0}.pdb amh-go.gro".format(proteinName))
    do("python2 ~/opt/Pdb2Gro.py crystal_structure.pdb amh-go.gro")
## ssweight
do("stride crystal_structure.pdb > ssweight.stride")
do("python2 ~/opt/script/stride2ssweight.py > ssweight")

## below used for old zim file
# do("grep -E 'CB|CA  GLY' crystal_structure.pdb > cbs.data")
# do("""awk '{if($9>15) print "1"; else if($9<-15) print "3"; else print "2"}'  cbs.data  > zim""")


# create "in" file
alpha_carbons = " ".join([str(i) for i in list(range(1, size*3+1, 3))])
beta_atoms = " ".join([str(i) for i in list(range(3, size*3+1, 3))])
oxygens = " ".join([str(i) for i in list(range(2, size*3+1, 3))])
last = " ".join(str(i) for i in [size*3-2, size*3-1, size*3])
if args.crystal:
    do("cp ~/opt/create_project_in_crystal_template.in {}_multi.in".format(proteinName))
else:
    if args.membrane:
        do("cp ~/opt/create_membrane_protein_folding_project_in_template.in {}_multi.in".format(proteinName))
    else:
        do("cp ~/opt/create_project_in_template.in {}_multi.in".format(proteinName))
with fileinput.FileInput("{}_multi.in".format(proteinName), inplace=True, backup='.bak') as file:
    for line in file:
        tmp = line.replace("ALPHA_CARBONS", alpha_carbons)
        tmp = tmp.replace("BETA_ATOMS", beta_atoms)
        tmp = tmp.replace("OXYGENS", oxygens)
        tmp = tmp.replace("PROTEIN", proteinName)
        tmp = tmp.replace("LAST", last)
        # tmp = BETA_ATOMS
        print(tmp, end='')

# print(alpha_carbons)
# create coord and data file

## start with crystal structure
do("python2 ~/opt/script/PDBToCoordinates.py crystal_structure crystal.coord")
do("python2 ~/opt/small_script/coord2data.py crystal.coord data.crystal -b")
if False:
    do("python2 ~/opt/script/PDBToCoordinates.py {0} {0}.coord".format(proteinName))
    do("python2 ~/opt/small_script/coord2data.py {0}.coord data.{0} -b".format(proteinName))
# if args.crystal:
#
# ## start with man made long string
# if not args.crystal:


# copy parameters
if args.globular:
    do("cp ~/opt/parameters/globular_parameters/* .")
if args.membrane:
    do("cp ~/opt/parameters/membrane/* .")
# task specific input

## frag memory generation
# if args.frag:
#     do("cp ~/opt/database/cullpdb_pc80_* .")
#     do("python2 ~/opt/script/prepFragsLAMW_index.py \
#     cullpdb_pc80_res3.0_R1.0_d160504_chains29712 %s.fasta 20 0" % proteinName)
# if args.frag:
#     do("cp ~/opt/database/cullpdb_pc80_* .")
#     do("python2 ~/opt/MultCha_prepFrags_index_HO.py \
#     cullpdb_pc80_res3.0_R1.0_d160504_chains29712 %s.fasta 20 0 95" % proteinName)
if args.frag:
    do("cp ~/opt/database/cullpdb_pc80_* .")
    do("python2 ~/opt/script/MultCha_prepFrags_index.py \
    cullpdb_pc80_res3.0_R1.0_d160504_chains29712 abeta42_2.fasta 20 1 9 > logfile")
