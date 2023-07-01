set terminal gif size 1000,1000 animate delay 5 optimize
set output 'animation.gif'

datafile = 'output.dat'

num_rows_per_block = 401

stats datafile nooutput
num_blocks = STATS_blocks

set xrange [-3000:3000]
set yrange [-3000:3000]
set zrange [-500:500]

set grid
set grid front lc rgb 'gray' lw 01
set grid xtics ytics ztics

set xyplane at -500

do for [row=1:num_rows_per_block] {
    start_row = row
    splot datafile every ::start_row::start_row using 2:3:4 \
                    with points pointtype 7 pointsize 0.5 notitle
}
