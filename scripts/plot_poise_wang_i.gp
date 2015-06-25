set datafile separator ","

set terminal pngcairo \
 transparent enhanced font "arial,14" fontscale 1.0 size 800, 600
set key center bottom box

set title ''
set ylabel 'u / u_max'
set xlabel 'y / H'
set output 'poise-wang-i.png'
plot \
    'data/poise/wang/implicit/0000/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.000', \
    'data/poise/wang/implicit/0001/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.001', \
    'data/poise/wang/implicit/0005/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.005', \
		'data/poise/wang/implicit/0010/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.010'
