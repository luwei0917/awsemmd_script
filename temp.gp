#!/usr/bin/env gnuplot

#set pointtype 5
#set xlabel "Annealing index"
#set ylabel "BestQ"

# set yrange [0.3:0.55]
n = 3
plot for [i = 0:n] '0_'.i w l t 'run '.i
pause - 1
