#!/bin/bash
if [ $# -ne 11 ]; then
	echo "Usage: exe PDB sim_label repeat T_sample Tmin Tmax dT binN k0 qmin qmax"
	exit
fi
pdbID_upper=$1
pdbID=$( echo $pdbID_upper | tr "[:upper:]" "[:lower:]" )
echo $pdbID

sim_label=$2
repeat=$3
T_sample=$4
T_min=$5
T_max=$6
dT=$7
binN_q=$8
k0=$9
L0=2
qmin=${10}
qmax=${11}

deltaq=$(echo "scale=4; $qmax - $qmin" | bc )
dq=$(echo "scale=4; $deltaq / ($repeat - 1) " | bc)
q=()
qlist=()
for (( j=0; j<$repeat; j++ )); do
    q0=$(echo "scale=4; $qmin + $dq * $j " | bc | awk '{printf "%f", $0}' )
    q=( ${q[*]} $q0 )
    if [ $j == 0 ]; then
	    qlist=( $qlist $q0 )
    else
	    qlist=( ${qlist[*]}, $q0 )
    fi
done
echo ${q[*]}
echo ${qlist[*]}


binN_E=100
tolerance=0.00001
n_T=$(echo "scale=0; ($T_max - $T_min)/$dT " | bc)
#n_T=$(( ($T_max - $T_min)/dT ))
echo $n_T

#q=( 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 )
#qlist=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

N=${#q[*]}
echo "# of simulation is $N "
#cp $qfile tmp
#qmin=${q[0]}
#qmax=${q[ $(($N-1)) ]}

cd $pdbID_upper

qfile=p_total
Ntotal=`cat $qfile | wc -l`
Nline=$(( $Ntotal / $N ))
echo "Nline=$Nline"

rm -rf simcount_file
for (( j=0; j<${#q[*]}; ++j )) ; do
	echo $Nline >> simcount_file
done

## build q_E_file
#awk '{ print $12 }' e_total > tmp
paste p_total e_total  > q_E_file

ss1=`~/opt/script/wham/Minmax $qfile`
qmax=`echo $ss1 | awk '{print $3}'`
qmin=`echo $ss1 | awk '{print $2}'`
echo "qmax: $qmax; qmin: $qmin"

ss=`~/opt/script/wham/Minmax e_total`
Emax=`echo $ss | awk '{print $3}'`
Emin=`echo $ss | awk '{print $2}'`
Eshift=`echo " ( $Emin + $Emax ) /2 " | bc`
echo $Emax, $Emin, $Eshift
sleep 2

for (( i_T=0 ; i_T<$n_T; ++i_T )) ; do
	T_target=$(( $T_min + ($T_max - $T_min)/$n_T*$i_T ))
	echo $T_target

	#generate input files
	sed "s/emin/$Emin/g" ~/opt/script/wham/Wham_generate_input.pl | sed "s/emax/$Emax/g" | sed "s/eshift/$Eshift/g" | \
	 sed "s/(xmin)/$qmin/" | sed "s/(xmax)/$qmax/" | sed "s/binN_q/$binN_q/" | sed "s/binN_E/$binN_E/" | sed "s/T_target/$T_target/" | sed "s/T_sample/$T_sample/" | \
	  sed "s/Nsample_list/$Nline/" | sed "s/qlist/${qlist[*]}/" | sed "s/k0/$k0/" | sed "s/L0/$L0/" | sed "s/q_E_file/q_E_file/" > tmp.pl
	chmod 755 tmp.pl
	./tmp.pl

	~/bin/wham -c counts_file -s simcount_file -u umbrella_potential_file -r $tolerance -x 4000 > whamout

	# convert Wham output to PMF files
	sed "s/emin/$Emin/g" ~/opt/script/wham/WhamPMF.pl | sed "s/emax/$Emax/g" | sed "s/eshift/$Eshift/g" | \
	 sed "s/(xmin)/$qmin/" | sed "s/(xmax)/$qmax/" | sed "s/binN_q/$binN_q/" | sed "s/binN_E/$binN_E/" | sed "s/T_target/$T_target/" | sed "s/T_sample/$T_sample/" | \
	  sed "s/sim_label/$sim_label/" > tmp.pl
	chmod 755 tmp.pl
	./tmp.pl
done
echo " $pdbID $sim_label $T_sample $T_min $T_max $dT $binN_q $k0"
