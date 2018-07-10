#!/usr/bin/env gnuplot
set term pdfcairo
set output 'q_ga.pdf'
#set pointtype 5
set xlabel "runs"
set ylabel "q value"
set title "GA"
set yrange [0:1]

plot "gagb.dat" u 1 w l
