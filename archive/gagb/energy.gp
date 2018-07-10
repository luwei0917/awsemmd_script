
set term pdfcairo
set output 'energy.pdf'
#set pointtype 5
set xlabel "million timesteps"
set ylabel "energy"

#set key outside

plot "wham.dat" u 1:2 w l, "wham_ga_rerun_with_gb_potential.dat" w l
#plot for [i=0:4] 'analysis/'.i.'/wham.dat' using ($1/10**6):2 w l smooth bezier t 'run '.i
#plot for [i=0:4] 'analysis/'.i.'/wham.dat' using ($1/10**6):($5<1000 ? $5 :1/0)  w l t 'run '.i
#smooth acsplines BEZIER  smooth {unique | frequency | cumulative | kdensity | csplines | acsplines | bezier | sbezier}
#plot for [i=0:9] 'wham.dat' using ($1/10**6):($5<1000 ? $5 :1/0)  w l t 'run '.i
