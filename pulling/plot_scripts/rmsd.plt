
set term pdfcairo
set output '../../results/rmsd_NUMBER.pdf'
#set pointtype 5
set xlabel "time"
set ylabel "rmsd"

plot "rmsd"  w l