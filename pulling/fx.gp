#!/usr/bin/env gnuplot
set term pdfcairo enhanced fontscale 0.5 size 5.00in, 3.00in
set output 'f_extension.pdf'
#set pointtype 5
set xlabel "Extension(Angstroms)"
set ylabel "Force(Kcal/mole-Angstrom)"
set title "2xov pulling"
#set yrange [-0.5:0.5]
set key outside

#plot "addforce.dat" u 2:3 w l
plot for [i=0:19] 'rerun_'.i.'/addforce.dat' u (($2)):(($3)) w l t 'run '.i \
,for [i=0:19] ''.i.'/addforce_back.dat' u 2:(-($3)) w l t 'run '.i
