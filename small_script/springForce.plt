#!/usr/bin/env gnuplot
set term pdfcairo
set output 'SpringForce.pdf'
#set pointtype 5
set xlabel "timesteps"
set ylabel "Spring Force"

#set yrange [-0.2:0.2]

plot "SpringForce.dat" u 1:5 w l
