set terminal pngcairo enhanced
set output 'min_energy_vs_h.png'
set title 'Minimum Energy vs Magnetic Field H'
set xlabel 'H (Magnetic Field)'
set ylabel 'minEnergy (Minimum Energy)'

# Definir os limites dos eixos
set xrange [0:5]
set yrange [-1:0]

# Definir os tics do eixo Y em intervalos de 0.1
set ytics 0.2
set xtics 0.5

# Gerar múltiplas linhas em um gráfico
plot for [i=1:8] sprintf('fort.%d', i) using 4:6 with lines title sprintf('Config %d', i)
