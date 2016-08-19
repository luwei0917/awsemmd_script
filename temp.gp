#!/usr/bin/env gnuplot

#set pointtype 5
#set xlabel "Annealing index"
#set ylabel "BestQ"

set key outside
# set yrange [0.3:0.55]
n = 3
plot for [j=0:n] for [i = 0:n] ''.j.'_'.i u 3 w lp pointtype (1+i)*(1+j) t 'run '.j.'\_'.i
pause - 1
