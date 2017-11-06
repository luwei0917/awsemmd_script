set term pdfcairo
set output 'list_of_max_q.pdf'
#set pointtype 5
set xlabel "Annealing index"
set ylabel "BestQ"

set yrange [0.3:0.55]

plot "/dev/stdin" u 1 w lp pt 7
