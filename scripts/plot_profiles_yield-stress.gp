set datafile separator ","

# set terminal pngcairo \
#   transparent enhanced font "arial,14" fontscale 1.0 size 1600, 1200
set key right box

set title 'Velocity profiles'
set ylabel 'u (lat unit / sec)'
set xlabel 'y (lat unit)'
# set output 'yield-stress.png'
plot \
    'data/bingham-mrt_5000_500_20_006_000008/prof-at-250_step-4000.dsv' \
			using 1:2 title 'nu=0.06, tau_y=8e-6', \
    'data/bingham-mrt_5000_500_20_006_000016/prof-at-250_step-4000.dsv' \
			using 1:2 title 'nu=0.06, tau_y=16e-6', \
		'data/bingham-mrt_5000_500_20_006_000032/prof-at-250_step-4000.dsv' \
			using 1:2 title 'nu=0.06, tau_y=32e-6', \
		'data/bingham-mrt_5000_500_20_006_000064/prof-at-250_step-4000.dsv' \
			using 1:2 title 'nu=0.06, tau_y=64e-6'
