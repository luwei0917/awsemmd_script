#!/usr/bin/env gnuplot
#set term pngcairo notransparent noenhanced fontscale 1 size 1080px, 30.00in
#set terminal epscairo
set term pngcairo notransparent
set output 'cf_extension.png'
#set pointtype 5
set ylabel "Extension(Angstroms)"
set xlabel "TimeSteps"
set title "2xov pulling"
#set yrange [-0.5:0.5]
set key outside

#plot "addforce.dat" u 2:3 w l

plot for [i=0:19] ''.i.'/addforce.dat' u 1:2 w l t 'run '.i
#plot for [i=0:19] ''.i.'/addforce.dat' u 2:3 w l t 'run '.i \
#,for [i=0:19] ''.i.'/addforce_back.dat' u 2:3 w l t 'run '.i
