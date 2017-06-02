
set term pdfcairo
set output 'results/cross_q.pdf'
set view map
#set size ratio .9

set key outside

#set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back
#set object 1 rect fc rgb "black" fillstyle solid 1.0

splot "results/cross_q" using 1:2:3 with points pointtype 5 pointsize 1 palette linewidth 20 t "crossq"
#splot "results/cross_q" using 1:2:3 palette linewidth 20 t "cross_q"
