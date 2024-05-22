set terminal png
set output 'plot01.png'
set title "Magnetização (m) vs Campo Magnético Externo (h)"
set xlabel "h"
set ylabel "m"
plot 'fort.1' using 1:2 with lines title 'N3', 'fort.2' using 1:2 with lines title 'N6', 'fort.3' using 1:2 with lines title 'N9', NaN with lines title 'Variância = 1', NaN with lines title 'Temperatura = 0.1'
