set datafile separator ","

# set terminal pngcairo \
#  transparent enhanced font "arial,14" fontscale 1.0 size 800, 600
set key center bottom box

set title ''
set ylabel 'u / u_max'
set xlabel 'x/W'
#set output 'yield-stress.png'

plot "data/poise_tauy-000000/ubar_profile.dsv" using 1:2 lt 1 lw 1 pt 1 title "tau = 0.00000", \
     "data/poise_tauy-000004/ubar_profile.dsv" using 1:2 lt 2 lw 1 pt 2 title "tau = 0.00004", \
     "data/poise_tauy-000008/ubar_profile.dsv" using 1:2 lt 3 lw 1 pt 3 title "tau = 0.00008", \
     "data/poise_tauy-000012/ubar_profile.dsv" using 1:2 lt 4 lw 1 pt 4 title "tau = 0.00012"
