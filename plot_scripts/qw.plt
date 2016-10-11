#!/usr/bin/env gnuplot
set term pdfcairo
set output 'qw.pdf'
#set pointtype 5
set xlabel "timesteps"
set ylabel "q value"

set yrange [0.2:0.9]

plot "wham.dat"  u 1:2 w l
