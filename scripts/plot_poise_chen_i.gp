set datafile separator ","

set terminal pngcairo \
 transparent enhanced font "arial,14" fontscale 1.0 size 800, 600
set key center bottom box

set title ''
set ylabel 'u / u_max'
set xlabel 'y / H'
set output 'poise-chen-i.png'
plot \
    'data/poise/chen/implicit/000000/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.00000', \
    'data/poise/chen/implicit/000004/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.00004', \
    'data/poise/chen/implicit/000008/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.00008', \
		'data/poise/chen/implicit/000012/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.00012'
