set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output '3d_output.png'

datafile = 'output.dat'

set title ''
set xlabel 'X axis [m]' font 'Verdana,24'
set ylabel 'Y axis [m]' font 'Verdana,24'
set zlabel 'Z axis [m]' font 'Verdana,24' rotate parallel

set xrange [*:*]
set yrange [*:*]
set zrange [*:*]

v_x = 0
v_y = 0
v_z = 1

set style line 1 lc rgb '#0000ff' lt 1 lw 3
set style increment user

set grid
set grid front lc rgb 'black' lw 2
set grid xtics ytics ztics

set xyplane at 0

set style arrow 1 size screen 0.02,15 filled lc rgb '#0000ff'

splot datafile using 2:3:4 index 0 with lines linestyle 1 notitle, \
      datafile using 2:3:4 index 0 with points linestyle 1 pointtype 7 pointsize 1 notitle
