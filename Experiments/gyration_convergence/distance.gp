set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output 'distance.png'

datafile = 'radius_error.dat'

set title ''
set xlabel 'log(N_{steps})' font ",28" offset 0.7
set ylabel 'log(|r_{fin} - r_0|)' font ",28" offset 0.8

set format x "10^{%L}"
set xrange [*:*]
set yrange [*:*]

set logscale y
set logscale x

set xtics font ",20"
set ytics font ",20"

set key font ",20"

plot datafile using 1:5 index 0 with linespoint pointtype 5 lc rgb "forest-green" pointsize 1.8 title "{/Symbol g} = 10" ,\
     datafile using 1:5 index 1 with linespoint pointtype 7 lc rgb "blue" pointsize 1.8 title "{/Symbol g} = 10^4" ,\
     datafile using 1:5 index 2 with linespoint pointtype 9 lc rgb "red" pointsize 1.8 title "{/Symbol g} = 10^6" 


