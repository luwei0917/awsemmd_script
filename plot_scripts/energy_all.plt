
set term pdfcairo
set output 'results/energy.pdf'
#set pointtype 5
set xlabel "2 million fs"
set ylabel "energy"

set key outside



plot for [i=0:9] 'analysis/'.i.'/wham.dat' using ($1/10**6):5 w l t 'run '.i
