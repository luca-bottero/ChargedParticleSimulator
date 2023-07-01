set terminal pngcairo size 1600,1200 enhanced font 'Verdana,16'
set output 'trajectories_2d.png'

datafile = 'output.dat'
b_field = "b_field.dat"

set title ''
set xlabel 'x axis'
set ylabel 'y axis'

plot b_field using 1:2:3:4 with vectors as 1 