set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output 'trajectories.png'

datafile = 'output.dat'

set title ''
set xlabel 'X axis [m]' rotate parallel
set ylabel 'Y axis [m]' rotate parallel
set zlabel 'Z axis [m]' rotate parallel

set xrange[-2000:2000]
set yrange[-2000:2000]

set xyplane at -30

splot for [i=0:*] datafile index i using 2:3:4 with linespoints pointtype 7 lc rgb "#FA0000F5" pointsize 0 notitle 
