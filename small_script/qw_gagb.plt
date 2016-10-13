#!/usr/bin/env gnuplot
set term pdfcairo
set output 'qw_gagb.pdf'
#set pointtype 5
set xlabel "(thousand)timesteps"
set ylabel "q value"

set yrange [0:0.6]

plot "q_ga.dat" w l,"q_gb.dat" w l
