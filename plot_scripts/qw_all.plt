#!/usr/bin/env gnuplot
set term pdfcairo
set output 'results/qw_all.pdf'
#set pointtype 5
set xlabel "million timesteps"
set ylabel "q value"

set key outside



#plot for [i=0:19] 'analysis/'.i.'/wham.dat' using ($1/10**6):2  w l smooth acsplines t 'run '.i
#smooth acsplines BEZIER  smooth {unique | frequency | cumulative | kdensity | csplines | acsplines | bezier | sbezier}
#plot for [i=0:9] 'wham.dat' using ($1/10**6):($5<1000 ? $5 :1/0)  w l t 'run '.i
plot for [i=0:19] 'analysis/'.i.'/wham.dat' using ($1/10**6):2  w l linetype i dashtype i smooth bezier t 'run '.i 
