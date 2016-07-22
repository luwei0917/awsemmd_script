

set term pdfcairo
set output out_file_name
#set pointtype 5
set xlabel "Annealing index"
set ylabel "BestQ"

set yrange [0.3:0.55]

plot in_file_name u 1 w lp pt 7
