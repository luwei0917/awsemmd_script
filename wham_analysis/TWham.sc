#!/bin/bash

if [ $# -ne 6 ]; then
	echo "Usage: exe PDB sim_label T_sample T_target binN qo_flag"
	exit
fi
pdbID_upper=$1
pdbID=$( echo $pdbID_upper | tr "[:upper:]" "[:lower:]" )
echo $pdbID

sim_label=$2
L0=4
k0=10000
T_sample=$3
T_target=$4
binN_q=$5
qo_flag=$6
binN_E=100
tolerance=0.00001

#q=(0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90)
#qlist=(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90)
qlist=( 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95)
#q=( 0.15 0.20 0.25 0.30 0.35 0.375 0.40 0.425 0.45 0.475 0.50 0.525 0.55 0.575 0.60 0.65 0.70 0.75 0.80 )
#qlist=(0.15, 0.20, 0.25, 0.30, 0.35, 0.375, 0.40, 0.425, 0.45, 0.475, 0.50, 0.525, 0.55, 0.575, 0.60, 0.65, 0.70, 0.75, 0.80 )

N=${#q[*]}
qmin=${q[0]}
qmax=${q[ $(($N-1)) ]}
cd $pdbID_upper

echo "# of simulations is $N."
Ntotal=`cat q_$sim_label | wc -l`
Nline=$(( $Ntotal / $N ))
echo "Nline=$Nline"

rm -rf simcount_file
for (( j=0; j<${#q[*]}; ++j )) ; do
	echo $Nline >> simcount_file
done

## build q_E_file
awk '{ print $12 }' energy_$sim_label > tmp
if [ $qo_flag -eq 1 ]; then
	paste qo_$sim_label tmp > q_E_file
else 
	paste q_$sim_label tmp > q_E_file
fi

ss=`Minmax tmp`
Emax=`echo $ss | awk '{print $3}'`
Emin=`echo $ss | awk '{print $2}'`
Eshift=`echo " ( $Emin + $Emax ) /2 " | bc`

echo $Emax, $Emin, $Eshift
#generate input files 
sed "s/emin/$Emin/g" ../Wham_generate_input.pl | sed "s/emax/$Emax/g" | sed "s/eshift/$Eshift/g" | \
	sed "s/(xmin)/$qmin/" | sed "s/(xmax)/$qmax/" | sed "s/binN_q/$binN_q/" | sed "s/binN_E/$binN_E/" | sed "s/T_target/$T_target/" | sed "s/T_sample/$T_sample/" | \
		sed "s/Nsample_list/$Nline/" | sed "s/qlist/${qlist[*]}/" | sed "s/k0/$k0/" | sed "s/L0/$L0/" | sed "s/q_E_file/q_E_file/" > tmp.pl
chmod 755 tmp.pl
./tmp.pl

wham -c counts_file -s simcount_file -u umbrella_potential_file -r $tolerance -x 4000 > whamout

# convert Wham output to PMF files
sed "s/emin/$Emin/g" ../WhamPMF.pl | sed "s/emax/$Emax/g" | sed "s/eshift/$Eshift/g" | \
	sed "s/(xmin)/$qmin/" | sed "s/(xmax)/$qmax/" | sed "s/binN_q/$binN_q/" | sed "s/binN_E/$binN_E/" | sed "s/T_target/$T_target/" | sed "s/T_sample/$T_sample/" | \
		sed "s/sim_label/$sim_label/" > tmp.pl
chmod 755 tmp.pl
./tmp.pl
