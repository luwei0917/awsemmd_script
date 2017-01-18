#!/bin/bash

# Usage ./run_frustratomer.sh PDBID #CHAINID
if [ $# -ne 2 ]; then
	echo "Usage ./run_frustratomer.sh PDBID #CHAINID"
	exit
fi

# Setup scripts and paths
#convpdbscript=/home/mc70/mmtsb/perl/convpdb.pl
#genereate input files
pdbcoords2lammpsscript=~/opt/script/PdbCoords2Lammps.sh
#parameters for running AWSEM simulation
parametersfolder=~/opt/para
#compile your own lammps exe
lmpfrustexec=~/bin/lmp_serial
strideexec=stride
stride2ssweightscript=~/opt/script/stride2ssweight.py

# make a directory
#mkdir -p frust$1
# go into the directory
#cp $1.pdb frust$1/
#cd frust$1
#if [ ! -f $1.pdb ];
#then
    # get the relevant pdb file
#    $pdbgetscript $1
    # clean the pdb and select the correct chain
#    $convpdbscript -chain $2 -model 1 -nohetero -renumber 1 -out generic -fixcoo $1.pdb > $1$2_cleaned.pdb
#    $convpdbscript -nohetero -renumber 1 -out generic -fixcoo $1.pdb > $1_cleaned.pdb
    # create the input files for AWSEM
#    $pdbcoords2lammpsscript $1_cleaned $1
    $pdbcoords2lammpsscript $1 $1
    # make the ssweight file
    $strideexec $1.pdb > ssweight.stride
    python2 $stride2ssweightscript > ssweight
#fi
# copy the necessary parameter files
cp $parametersfolder/* .
if [ -e ./charge_on_residues.dat ]; then
	cp ./charge_on_residues.dat ./ -v
else
	echo "0" > charge_on_residues.dat
	echo "No charge files, generated an empty one!"
	sleep 2
fi

sed -i "s/run\t\t10000/run 0/" $1.in
# run the calculation
$lmpfrustexec -in $1.in
