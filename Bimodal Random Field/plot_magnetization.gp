set terminal pngcairo enhanced
set output 'magnetization_vs_h.png'
set title 'Magnetization M(i)/N vs Magnetic Field H'
set xlabel 'H (Magnetic Field)'
set ylabel 'M(i)/N (Magnetization per Spin)'

# Gerar múltiplas linhas em um gráfico
plot for [i=1:8] sprintf('fort.%d', i) using 4:5 with lines title sprintf('Config %d', i)
