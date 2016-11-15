#!/usr/bin/env gnuplot
set term pdfcairo
set output 'compare.pdf'
#set pointtype 5
set xlabel "Qw based on ga"
set ylabel "Free energy(kcal/mol)"
set title "GAGB free energy with ga qw"
#set yrange [-0.5:0.5]

plot "ga_seq.dat" w l, "gb_seq.dat" w l
