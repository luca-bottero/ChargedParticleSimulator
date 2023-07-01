set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output 'dipole.png'

datafile = 'output.dat'

set title ''
set xlabel 'X axis [R_{Earth}}]' rotate parallel
set ylabel 'Y axis [R_{Earth}}]' rotate parallel
set zlabel 'Z axis [R_{Earth}}]' rotate parallel

set xrange [-3:3]
set yrange [-3:3]
set zrange [-3:3]

set xyplane at -1.5

set view equal

set parametric
set isosamples 50,50

set urange [-pi/2:pi/2]
set vrange [0:2*pi]
set ztics nomirror -1.0,0.25,1.0

splot for [i=0:*] datafile index i using 2:3:4 with linespoints pointtype 7 pointsize 0 notitle, \
        cos(u)*cos(v), cos(u)*sin(v), sin(u) with lines lc rgb "#F0000000" notitle
