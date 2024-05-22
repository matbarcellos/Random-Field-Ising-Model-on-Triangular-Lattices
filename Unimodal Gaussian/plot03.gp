set terminal png
set output 'plot03.png'
set title "Magnetização (m) vs Campo Magnético Externo (h)"
set xlabel "h"
set ylabel "m"
plot 'fort.24' using 1:2 with lines title 'N3', 'fort.25' using 1:2 with lines title 'N6', 'fort.26' using 1:2 with lines title 'N9', NaN with lines title 'Variância = 0.5', NaN with lines title 'Temperatura = 0.1'
