#!/bin/bash
if [ $# -ne 13 ]; then
	echo "Usage: exe PDB 'sim_list' repeat 'T_sample_list' Tmin Tmax dT binN k0 qmin qmax new_sim_label qname(q,qo,qdiff)"
	exit
fi
pdbID_upper=$1
pdbID=$( echo $pdbID_upper | tr "[:upper:]" "[:lower:]" )
echo $pdbID
sim_label_list=( $2 )
repeat=$3
Tsample_list=( $4 )
T_min=$5
T_max=$6
dT=$7
binN_q=$8
k0=$9
qmin=${10}
qmax=${11}
new_sim_label=${12}
qname=${13}
echo $qname

if [ ${#sim_label_list[*]} -ne ${#Tsample_list[*]} ]; then
	echo "Warning: sim list has different length as T-sample list! Fix it first."
	echo "Exiting..."
	exit
fi
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
n_T=$(( ($T_max - $T_min)/dT ))
echo $n_T

L0=2

N=${#q[*]}
qmin=${q[0]}
qmax=${q[ $(($N-1)) ]}
cd $pdbID_upper

#generate simcount_file_mult
echo "# of simulations is $N."
rm -f simcount_file_mult E_tmp q_tmp
rm -f $qname\_$new_sim_label energy_$new_sim_label
for (( i=0; i<${#sim_label_list[*]}; ++i )); do
	sim_label=${sim_label_list[$i]}
	qfile=$qname\_$sim_label
	Ntotal=`cat $qfile | wc -l`
	Nline=$(( $Ntotal / $N ))
	echo "Nline=$Nline"

	for (( j=0; j<${#q[*]}; ++j )) ; do
		echo $Nline >> simcount_file_mult
	done
	## build q_E_file
	awk '{ print $NF }' energy_$sim_label >> E_tmp
	cat $qfile >> q_tmp
	#fi

	cat $qfile      >> $qname\_$new_sim_label
	cat energy_$sim_label >> energy_$new_sim_label
	#cat all_contact_$sim_label >> all_contact_$new_sim_label
done

## build q_E_file
paste q_tmp E_tmp > q_E_file
newqfile=$qname\_$new_sim_label

ss1=`~/opt/script/wham/Minmax $newqfile`
qmax=`echo $ss1 | awk '{print $3}'`
qmin=`echo $ss1 | awk '{print $2}'`
echo "qmax: $qmax; qmin: $qmin"

ss=`~/opt/script/wham/Minmax E_tmp`
Emax=`echo $ss | awk '{print $3}'`
Emin=`echo $ss | awk '{print $2}'`
Eshift=`echo " ( $Emin + $Emax ) /2 " | bc`
echo $Emax, $Emin, $Eshift
sleep 4

for (( i_T=0 ; i_T<$n_T; ++i_T )) ; do
	T_target=$(( $T_min + ($T_max - $T_min)/$n_T*$i_T ))
	echo $T_target

	#generate input files
	rm -f tmp_umbrella
	for (( i=0; i<${#Tsample_list[*]}; ++i )); do
		T_sample=${Tsample_list[$i]}
		sed "s/emin/$Emin/g" ~/opt/script/wham/Wham_generate_input.pl | sed "s/emax/$Emax/g" | sed "s/eshift/$Eshift/g" | \
		 sed "s/(xmin)/$qmin/" | sed "s/(xmax)/$qmax/" | sed "s/binN_q/$binN_q/" | sed "s/binN_E/$binN_E/" | sed "s/T_target/$T_target/" | sed "s/T_sample/$T_sample/" | \
		  sed "s/Nsample_list/$Nline/" | sed "s/qlist/${qlist[*]}/" | sed "s/k0/$k0/" | sed "s/L0/$L0/" | sed "s/q_E_file/q_E_file/" > tmp.pl
		chmod 755 tmp.pl
		./tmp.pl
		cat umbrella_potential_file >> tmp_umbrella
	done
	~/bin/wham -c counts_file -s simcount_file_mult -u tmp_umbrella -r $tolerance -x 4000 > whamout

	# convert Wham output to PMF files
	sed "s/emin/$Emin/g" ~/opt/script/wham/WhamPMF.pl | sed "s/emax/$Emax/g" | sed "s/eshift/$Eshift/g" | \
	 sed "s/(xmin)/$qmin/" | sed "s/(xmax)/$qmax/" | sed "s/binN_q/$binN_q/" | sed "s/binN_E/$binN_E/" | sed "s/T_target/$T_target/" | \
	  sed "s/sim_label/$new_sim_label/" > tmp.pl
	chmod 755 tmp.pl
	./tmp.pl
done
echo " $pdbID $new_sim_label $T_target $T_min $T_max $dT $binN_q"
