set terminal gif size 1600,1200 animate delay 2 optimize
set output '3d_animation.gif'

datafile = 'output.dat'

num_rows_per_block = 20001

stats datafile nooutput
num_blocks = STATS_blocks

set xrange [-3000:3000]
set yrange [-3000:3000]
set zrange [-5:400]

set grid
set grid front lc rgb 'gray' lw 01
set grid xtics ytics ztics

set xyplane at -5

do for [row=1:num_rows_per_block] {
    start_row = row
    splot datafile every ::start_row::start_row using 2:3:4 \
                    with points pointtype 7 lt rgb "blue" pointsize 0.5 notitle
}
