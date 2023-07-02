set terminal gif size 1000,1000 animate delay 0.01 optimize
set output 'animation.gif'

datafile = 'output.dat'

num_rows_per_block = 1000

stats datafile nooutput
num_blocks = STATS_blocks

set title ''
set xlabel 'X axis [R_{Earth}]' rotate parallel
set ylabel 'Y axis [R_{Earth}]' rotate parallel
set zlabel 'Z axis [R_{Earth}]' rotate parallel

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

do for [row=1:num_rows_per_block] {
    start_row = row
    splot datafile every ::0::200*start_row using 2:3:4 \
                    with points pointtype 7 pointsize 0.5 notitle,\
                    cos(u)*cos(v), cos(u)*sin(v), sin(u) with lines lc rgb "#F0000000" notitle

}
