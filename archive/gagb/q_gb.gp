#!/usr/bin/env gnuplot
set term pdfcairo
set output 'q_gb.pdf'
#set pointtype 5
set xlabel "runs"
set ylabel "q value"
set title "GB"
set yrange [0:1]

plot "gagb.dat" u 2 w l
