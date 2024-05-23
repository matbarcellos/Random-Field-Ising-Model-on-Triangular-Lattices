# plot_all.gp
set terminal pngcairo enhanced size 800,600
set output 'combined_plot.png'

set title "M vs H/|J|"
set xlabel "H/|J|"
set ylabel "M"

set grid
set key outside

plot 'H_0_1.dat' using 1:2 with lines title 'Dataset 1', \
     'H_0_2.dat' using 1:2 with lines title 'Dataset 2', \
     'H_0_3.dat' using 1:2 with lines title 'Dataset 3', \
     'H_0_4.dat' using 1:2 with lines title 'Dataset 4', \
     'H_0_5.dat' using 1:2 with lines title 'Dataset 5'