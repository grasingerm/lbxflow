set datafile separator ","

set key center bottom box

set title ''
set ylabel 'u / u_max'
set xlabel 'y / H'
set output 'poise-wang-e.png'
plot \
    'data/poise/wang/explicit/0001/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.001', \
    'data/poise/wang/explicit/0005/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.005', \
		'data/poise/wang/explicit/0010/ux_profile.dsv' \
      using 1:2 title 'yield stress = 0.010'
