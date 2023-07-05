set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output 'phase.png'

datafile = 'radius_error.dat'

set title ''
set xlabel 'log(dt/dt_0)' font ",28"
set ylabel 'log({/Symbol f})' font ",28" offset 0.8

set xrange [*:*]
set yrange [*:*]

set logscale y
set logscale x

set key font ",20"

plot datafile using 2:4 index 0 with linespoint pointtype 5 lc rgb "forest-green" pointsize 1.5 title "{/Symbol g} = 10" ,\
     datafile using 2:4 index 1 with linespoint pointtype 7 lc rgb "blue" pointsize 1.5 title "{/Symbol g} = 10^4" 

