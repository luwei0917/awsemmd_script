#!/bin/bash

export s=$HOME/opt

if [[ ! $# -eq 2 ]]
then
        echo
        echo "> $0 pdb_id project_name"
        echo
        exit
fi

pdb_file=$1
output_file=$2

echo $pdb_file
echo $output_file

python2.7 $s/script/PDBToCoordinates.py $pdb_file $output_file".coord"
python2.7 $s/small_script/coord2data.py $output_file".coord" "data."$output_file -b
