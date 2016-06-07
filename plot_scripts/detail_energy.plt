set term pdfcairo
set output 'detail_energy.pdf'
#set pointtype 5
set xlabel "time"
set ylabel "energy"

#plot "energy.dat"  u 1:9 w l,"energy.dat" u 1:10  w l ,"energy.dat" u 1:11  w l
set key outside
plot for [col=2:20] 'energy.dat' using 0:col with lines title columnheader
