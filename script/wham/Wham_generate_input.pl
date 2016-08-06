#!/usr/bin/perl 
# Weihua Jun.07
# 2dpotential, modified from perlgsWham.
# create files for Wham. counts_file/ umbrella_potential, simcount_file is created manually.
# target temperature is determined by what you used to generatre umbrella_potential at about line 57

$kb=0.001987 ;
$Emin= emin - ( eshift ) ; $Emax = emax - ( eshift ) ; #17.703240;
$xmin = (xmin) ; $xmax = (xmax) ; 
$bin_n_x = binN_q ; $bin_n_E = binN_E ;
$beta0  = 1.0/($kb*T_target) ;
$betaj  = 1.0/($kb*T_sample) ;
$k_bias = k0 ;
$L_bias = L0 ;
#@beta = ( 1.0/T_sample ) ;
@N = ( Nsample_list ); #number of states in each walker file
@q = ( qlist ) ;

$size = @q ;
#$size = @beta ;
$xbin_size =($xmax-$xmin)/$bin_n_x;
$Ebin_size =($Emax-$Emin)/$bin_n_E;

foreach $j(0..$bin_n_x*$bin_n_E-1){$count[$j]=0;}
open IN, "q_E_file" or die "can't open data file: $!"; 
foreach (<IN>){
	@line =split(' ');
	$x =$line[0]; $E =$line[1] - ( eshift ) ; 

        if ($x<=$xmin){ $k1 = 0; }
        elsif ($x>=$xmax){ $k1 = $bin_n_x -1 ; }
	else { $k1 = int(($x - $xmin) / $xbin_size ) ; }
        if ($E<=$Emin){ $k2 = 0; }
        elsif ($E>=$Emax){ $k2 = $bin_n_E -1 ; }
	else { $k2 = int(($E - $Emin) / $Ebin_size ) ; }
	#print "$x $E $k1 $k2\n";

	$k = $k1 + $k2*$bin_n_x;
	$count[$k]++;
}
close IN;

####################################################################
#generate wham counts file
$i=0;
open OUT, ">./counts_file";
foreach $j (0..$bin_n_x*$bin_n_E-1){
	print OUT "$count[$j]\n";
	$i++;
}
print "counts file has $i line\n";
close OUT;

##################################################################
#create umbrella_potential file
open OUT2, ">./umbrella_potential_file";
foreach $j (0..$size-1){
	open IN, "counts_file" or die "can't open data file: $!"; 
	$c_max=0;
	$i=0;
	$qj = $q[$j] ;
	#print "$qj ", "$betaj ", "$beta0\n" ;
	#print "$k_bias ", "$L_bias\n" ;
	foreach (<IN>){
		$x_index = ($i)%$bin_n_x;
		$E_index = ($i-$x_index)/$bin_n_x;
		$E = $Emin + $Ebin_size * ($E_index + 0.5);
		$q = $xmin + $xbin_size * ($x_index + 0.5);
		$V = $k_bias * ($q - $qj)**$L_bias ;
		#print "$q $E $V\n" ;
		$c_i = exp(-($betaj-$beta0) * $E ) * exp( -$betaj * $V ) ;

		if($c_max<$c_i){$c_max =$c_i;}
		print OUT2 "$c_i ";
		$i++;
	}
	print OUT2 "\n";
	close IN;
	print "cij_max: $c_max\n";
}
close OUT2;
