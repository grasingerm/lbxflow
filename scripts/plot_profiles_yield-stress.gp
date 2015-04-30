set datafile separator ","

set terminal pngcairo \
  transparent enhanced font "arial,14" fontscale 1.0 size 800, 600
set key center bottom box

set title ''
set ylabel 'velocity (lat unit / sec)'
set xlabel 'distance across channel (lat unit)'
set output 'yield-stress.png'
plot \
    'data/bingham-mrt_5000_500_20_006_000008/prof-at-250_step-5000.dsv' \
      using 1:2 title 'yield stress = 0.000', \
    'data/bingham-mrt_5000_500_20_006_004000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'yield stress = 0.004', \
    'data/bingham-mrt_5000_500_20_006_008000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'yield stress = 0.008', \
		'data/bingham-mrt_5000_500_20_006_016000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'yield stress = 0.016', \
		'data/bingham-mrt_5000_500_20_006_032000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'yield stress = 0.032', \
		'data/bingham-mrt_5000_500_20_006_064000/prof-at-250_step-5000.dsv' \
			using 1:2 title 'yield stress = 0.064'
