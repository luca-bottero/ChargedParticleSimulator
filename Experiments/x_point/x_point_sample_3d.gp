set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output 'x_point_sample_3d.png'

datafile = 'output.dat'

set title ''
set xlabel 'X axis [m]' rotate parallel
set ylabel 'Y axis [m]' rotate parallel
set zlabel 'Z axis [m]' rotate parallel

set xrange[-1000:1000]
set yrange[-1000:1000]

set xyplane at -30

splot for [i=0:10] datafile index i using 2:3:4 with linespoints pointtype 7 pointsize 0 notitle 
