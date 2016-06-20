
set term pdfcairo
set output 'results/energy.pdf'
#set pointtype 5
set xlabel "million timesteps"
set ylabel "energy"

set key outside



plot for [i=0:9] 'analysis/'.i.'/wham.dat' using ($1/10**6):($5<1000 ? $5 :1/0)  w l t 'run '.i
#smooth acsplines BEZIER  smooth {unique | frequency | cumulative | kdensity | csplines | acsplines | bezier | sbezier}
