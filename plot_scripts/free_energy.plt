
set term pdfcairo
set output 'free_energy.pdf'
#set pointtype 5
set xlabel "q"
set ylabel "F(kcal/mol)"
if (!exists("number_of_run")) number_of_run=19
if (!exists("temperature")) temperature=300

n = 1
dt = 100
#plot "pmf-".temperature.".dat"  u 2:3 w l
#plot for [i=0:3] "pmf-".(400+i*50)."-ga.dat"  u 2:3 w l dt 2 , "pmf-450-ga.dat"  u 2:3 w l dt 2 , "pmf-500-ga.dat"  u 2:3 w l dt 2 ,"pmf-400.dat"  u 2:3 w l, "pmf-450.dat"  u 2:3 w l, "pmf-500.dat"  u 2:3 w l
#plot for [i=0:9] 'analysis/'.i.'/wham.dat'
plot for [i=0:n] "ga/pmf-".(temperature+i*dt).".dat"  u 2:3 w l lt (i+1) dt 1 t "ga T ".(temperature+i*dt) ,\
for [i=0:n] "gb/pmf-".(temperature+i*dt).".dat"  u 2:3 w l lt (i+1) dt 2 t "gb T ".(temperature+i*dt)
