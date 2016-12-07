#!/usr/bin/env gnuplot
set term pdfcairo enhanced notransparent fontscale 0.5 size 5.00in, 3.00in
set output 'oneFreeEnergy.pdf'
#set pointtype 5
set xlabel "Extension(Angstroms)"
set ylabel "Free Energy(Kcal/mole)"
set title "2xov pulling"
#set yrange [-0.5:0.5]
set key outside

#plot "addforce.dat" u 2:3 w l
plot 'pmf-300.dat' u 2:3 w l 
