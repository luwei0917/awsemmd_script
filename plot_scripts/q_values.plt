set term pdfcairo
set output 'q_values.pdf'
#set pointtype 5
set xlabel "Annealing index"
set ylabel "BestQ"

plot "ha_q" w lp pt 7,"he_q" w lp pt 5
