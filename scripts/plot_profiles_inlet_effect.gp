set datafile separator ","

set terminal pngcairo \
    transparent enhanced font "arial,14" fontscale 1.0 size 1600, 1200
set key right box

set title 'Velocity profiles'
set ylabel 'u (lat unit / sec)'
set xlabel 'y (lat unit)'
set output 'inlet-effect.png'
plot \
    'data/prof-at-100_step-10000.dsv' using 1:2 title 'x=100', \
    'data/prof-at-250_step-10000.dsv' using 1:2 title 'x=250', \
    'data/prof-at-500_step-10000.dsv' using 1:2 title 'x=500'
