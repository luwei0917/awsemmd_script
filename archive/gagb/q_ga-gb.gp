#!/usr/bin/env gnuplot
set term pdfcairo
set output 'q_ga-gb.pdf'
#set pointtype 5
set xlabel "runs"
set ylabel "q value"
set title "GA-GB"
set yrange [-0.5:0.5]

plot "gagb.dat" u ($1-$2) w l
