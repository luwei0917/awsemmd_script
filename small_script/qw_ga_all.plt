#!/usr/bin/env gnuplot
set term pdfcairo
set output 'results/qw_ga_all.pdf'
#set pointtype 5
set xlabel "thousand timesteps"
set ylabel "q value"

set key outside
if (!exists("number_of_run")) number_of_run=19


#plot for [i=0:number_of_run] 'analysis/'.i.'/wham.dat' using ($1/10**6):2  w l smooth acsplines t 'run '.i
#smooth acsplines BEZIER  smooth {unique | frequency | cumulative | kdensity | csplines | acsplines | bezier | sbezier}
#plot for [i=0:9] 'wham.dat' using ($1/10**6):($5<1000 ? $5 :1/0)  w l t 'run '.i
plot for [i=0:number_of_run] 'analysis/'.i.'/q_ga.dat' w l linetype i dashtype i smooth bezier t 'run '.i
