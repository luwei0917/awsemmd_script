#!/usr/bin/perl 
# Weihua Jun.07
# 2dpotential
# using wham output (unbiased distribution) to generate PMF.

$kb=0.001987 ;
$Emin= -517.62451114169834909 - ( -431 ) ; $Emax = -344.90577690292542457 - ( -431 ) ; #17.703240;
$xmin = 0.107 ; $xmax = 0.763 ; 
$bin_n_x = 50 ; $bin_n_E = 100 ;
$beta0  = 1.0/($kb*365) ;
$betaj  = 1.0/($kb*350) ;

$xbin_size =($xmax-$xmin)/$bin_n_x;
$Ebin_size =($Emax-$Emin)/$bin_n_E;

open OUT, ">./T089_365_pmf.dat";

open IN, "./whamout" or die "can't open data file: $!"; 
foreach (<IN>){
	@line =split(' ');
	push(@prob, $line[0]);
}
$size = @prob;
if ($size != $bin_n_x*$bin_n_E){print "Error: Size of whamout mismatches with total bin number.\n";}
close IN;

foreach $i (0..$bin_n_x-1){$p_x[$i]=0;}
$p_check=0 ; 
foreach $i (0..$size-1){
	$x_index = ($i)%$bin_n_x;
	$p_x[$x_index] += $prob[$i]/$xbin_size; #transform probability to density
	$p_check+= $prob[$i]; 
}
print "check $p_check\n";
print "beta=$beta0\n";
#$p_check=0 ; 
#foreach $i (0..$bin_n_x-1){$p_check+=$p_x[$i]*$xbin_size;}
#print "check again $p_check\n";

foreach $i (0..$bin_n_x-1){
	$x=$xmin+($i+0.5)*$xbin_size;
	if ($p_x[$i] != 0){ 
		$pmf = -log($p_x[$i])/$beta0 ; 
		print OUT "$x $pmf\n";
	}
	else { 
		print OUT "$x inf\n";
	}
}
close OUT;
