set autoscale
set xtic auto
set ytic auto
set xlabel "Angle of attack (rad)"
set ylabel "C_mq"

plot "all_aero_coefs.dat" using ($2 * 3.14 / 180):20 with linespoints