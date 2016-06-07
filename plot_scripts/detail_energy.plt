set term pdfcairo
set output 'detail_energy.pdf'
#set pointtype 5
set xlabel "time"
set ylabel "energy"

plot "energy.dat"  u 1:5 w l,"energy.dat" u 1:6  w l 
