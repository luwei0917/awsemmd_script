#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import myPersonalFunctions
import fileinput
from small_script.myFunctions import *

if(platform.system() == 'Darwin'):  # Mac system (local machine)
    OPENAWSEM_LOCATION = "/Users/weilu/openmmawsem/"
elif(platform.system() == 'Linux'):
    OPENAWSEM_LOCATION = '/projects/pw8/wl45/openmmawsem/'
else:
    print("system unknown")
sys.path.insert(0, OPENAWSEM_LOCATION)
from openmmawsem import *


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
parser.add_argument("--hybrid", action="store_true", default=False)


args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

proteinName = pdb_id = args.protein
# chain='A'
chain='ABC'
pdb = f"{pdb_id}.pdb"

# print(args)
with open('commandline_args.txt', 'w') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


# for compute Q
input_pdb_filename, cleaned_pdb_filename = prepare_pdb("crystal_structure.pdb", chain)
ensure_atom_order(input_pdb_filename)
# get fasta, pdb, seq file ready
getSeqFromCleanPdb(input_pdb_filename, chains=chain)
if not args.crystal:
    do("~/opt/fasta2pdb.py "+proteinName)
else:
    do(f"cp crystal_structure.pdb {pdb}")

input_pdb_filename, cleaned_pdb_filename = prepare_pdb(pdb, chain)
ensure_atom_order(input_pdb_filename)

if args.globular:
    os.system("cp /Users/weilu/opt/parameters/globular_parameters/burial_gamma.dat .")
    os.system("cp /Users/weilu/opt/parameters/globular_parameters/gamma.dat .")

do("python2 ~/opt/Pdb2Gro.py crystal_structure.pdb amh-go.gro")

## ssweight
do("stride crystal_structure.pdb > ssweight.stride")
do("python2 ~/opt/script/stride2ssweight.py > ssweight")

# below used for zim and zimPosition file
if args.membrane or args.hybrid:
    do("grep -E 'CB|CA  GLY' crystal_structure-cleanded.pdb > cbs.data")
    do("""awk '{if($9>15) print "1"; else if($9<-15) print "3"; else print "2"}'  cbs.data  > zimPosition""")
    if args.crystal:
        create_zim(f"crystal.seq")


if args.frag:
    do("cp ~/opt/database/cullpdb_pc80_* .")
    do("python2 ~/opt/script/MultCha_prepFrags_index.py \
    cullpdb_pc80_res3.0_R1.0_d160504_chains29712 %s.fasta 20 1 9 > logfile" % proteinName)
