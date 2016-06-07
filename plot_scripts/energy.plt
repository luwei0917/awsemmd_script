
set term pdfcairo
set output '../../results/energy_NUMBER.pdf'
#set pointtype 5
set xlabel "time"
set ylabel "energy"

plot "wham.dat"  u 1:5 w l
