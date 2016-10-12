#!/usr/bin/env gnuplot
set term pdfcairo
set output 'qw.pdf'
#set pointtype 5
set xlabel "timesteps"
set ylabel "q value"

set yrange [0.2:0.9]

plot "q_ga.dat" w l,"q_gb.dat" w l
