set autoscale
set xtic auto
set ytic auto
set xlabel "Angle of attack (rad)"
set ylabel "Cm,y about cg"

f(x) = m*x+b
m=-1.0; b=-0.01
fit f(x) "all_aero_coefs.dat" using (3.1415/180*$5):15 via m, b

plot "all_aero_coefs.dat" using (3.1415/180*$5):15:($15 / 10) with yerrorbars, m*x+b title "Linear fit"