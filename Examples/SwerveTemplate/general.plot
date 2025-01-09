set autoscale
set xtic auto
set ytic auto
set xlabel "Center of gravity / total length"
set ylabel "C_mq"

plot "all_aero_coefs.dat" using ($2 / 2.75):20 with linespoints