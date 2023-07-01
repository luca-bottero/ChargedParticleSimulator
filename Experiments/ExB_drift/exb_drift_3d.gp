set terminal pngcairo size 1600,1200 enhanced font 'Verdana,22'
set output 'drift_trajectories_3d.png'

datafile = 'output.dat'

set title ''
set xlabel 'X axis [m]' font 'Verdana,34' rotate parallel
set ylabel 'Y axis [m]' font 'Verdana,34' rotate parallel
set zlabel 'Z axis [m]' font 'Verdana,34' rotate parallel offset -1


set xrange [*:*]
set yrange [*:*]
set xyplane at 0

set label "E_{0,s}" right at first -0.01, second -0.0 font 'Verdana,44'
set label "E_{0,e}" right at first -0.007, second -0.0055 font 'Verdana,44'
set label "E_{0,w}" right at first -0.01, second 0.008 font 'Verdana,44'

splot datafile index 2 using 2:3:4 with points pointtype 7 lt rgb "blue" pointsize 1 notitle,\
      datafile index 1 using 2:3:4 with points pointtype 7 lt rgb "green" pointsize 1 notitle,\
      datafile index 0 using 2:3:4 with points pointtype 7 lt rgb "red" pointsize 1 notitle 

