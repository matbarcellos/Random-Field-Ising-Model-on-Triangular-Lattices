set terminal png
set output 'plot02.png'
set title "Magnetização (m) vs Campo Magnético Externo (h)"
set xlabel "h"
set ylabel "m"
plot 'fort.14' using 1:2 with lines title 'N3', 'fort.15' using 1:2 with lines title 'N6', 'fort.16' using 1:2 with lines title 'N9', NaN with lines title 'Variância = 1.5', NaN with lines title 'Temperatura = 0.1'
