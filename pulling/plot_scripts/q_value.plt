
set term pdfcairo
set output '../../results/q_NUMBER.pdf'
#set pointtype 5
set xlabel "time"
set ylabel "q_value"

plot "wham.dat"  u 1:2 w l