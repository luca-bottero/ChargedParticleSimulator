set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output 'norm_trajectories_plot.png'

datafile = 'normalized_trajectories.dat'

set title ''
set xlabel 'X axis [a.u.]' font ",36"
set ylabel 'Y axis [a.u.]' font ",36"

set xrange [*:*]
set yrange [*:*]

set size square
set grid

plot for[i=0:*] datafile index i using 2:3 with points pointtype 7 pointsize 1 notitle

