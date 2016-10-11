
file1 = in_file_name . "_ha.dat"
file2 = in_file_name . "_he.dat"
file3 = in_file_name . "_lp.dat"
file4 = in_file_name . "_he_lp.dat"
set term pdfcairo
set output out_file_name
#set pointtype 5
set xlabel "Annealing index"
set ylabel "BestQ"

#set yrange [0.3:0.55]
set yrange [0.2:1]
plot file1 u 1 w lp pt 7 t "Homolog allowed",\
 file2 u 1 w lp pt 7 t "Homolog excluded", \
 file3 u 1 w lp pt 7 t "Lowest Potential second run", \
 file4 u 1 w lp pt 7 t "LP from HE"
