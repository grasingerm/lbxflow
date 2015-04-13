set datafile separator ","

set terminal pngcairo \
  transparent enhanced font "arial,14" fontscale 1.0 size 1600, 1200
set key right bottom box

set title 'Velocity profiles'
set ylabel 'u (lat unit / sec)'
set xlabel 'y (lat unit)'
# set output 'yield-stress.png'
plot \
    'data/bingham-mrt_5000_500_20_006_004000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'tau_y=4e-3', \
    'data/bingham-mrt_5000_500_20_006_008000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'tau_y=8e-3', \
		'data/bingham-mrt_5000_500_20_006_016000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'tau_y=16e-3', \
		'data/bingham-mrt_5000_500_20_006_032000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'tau_y=32e-3', \
		'data/bingham-mrt_5000_500_20_006_064000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'tau_y=64e-3'
