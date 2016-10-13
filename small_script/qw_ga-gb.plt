#!/usr/bin/env gnuplot
set term pdfcairo
set output 'qw_ga-gb.pdf'
#set pointtype 5
set xlabel "(thousand)timesteps"
set ylabel "q value"

set yrange [-0.2:0.2]

plot "q_gagb.dat" u ($1-$2) w l
