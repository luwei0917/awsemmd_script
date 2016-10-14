#!/usr/bin/env gnuplot
set term pdfcairo
set output 'springForce_smooth.pdf'
#set pointtype 5
set xlabel "timesteps"
set ylabel "Spring Force"

#set yrange [-0.2:0.2]

plot "springForce.dat" u 1:5 w l smooth bezier
