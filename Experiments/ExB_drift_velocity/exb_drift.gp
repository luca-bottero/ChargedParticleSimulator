set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output 'plot.png'

datafile = 'output.dat'

set title ''
set xlabel 'x axis'
set ylabel 'y axis'

set xrange [*:*]
set yrange [*:*]

set logscale xy

plot datafile using 1:2 index 0 with points pointtype 7 pointsize 1 notitle

