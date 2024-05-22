set terminal pngcairo enhanced
set output 'grafico.png'
set xlabel 'Campo Magnético Externo'
set ylabel 'Magnetização'
set title 'Comparação de Magnetização com e sem Campos Aleatórios'
set grid
set key bottom right  # Posiciona a legenda na parte inferior direita

plot 'fort.1' using 1:2 with lines title 'Sem Campos Aleatórios', \
     'fort.15' using 1:2 with lines title 'Com Campos Aleatórios', \
     NaN title 'N = 15'

pause -1
