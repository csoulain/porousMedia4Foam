reset
set terminal postscript enhanced
set size square 0.65,0.65
set output "Porosity.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Porosity"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65



plot "postProcessing/graphCell/3.149928e+08/line.xy" using 1:2 with lines lw 5 title "10 Y", \
     "postProcessing/graphCell/3.1499928e+09/line.xy" using 1:2 with lines lw 5 title "100 Y", \
     "postProcessing/graphCell/3.78e+09/line.xy" using 1:2 with lines lw 5 title "120 Y", \
     "plots/MIN3P_Data_Paper/10Y.txt"  using 1:3	title "MIN3P 10Y", \
     "plots/MIN3P_Data_Paper/100Y.txt" using 1:3	title "MIN3P 100Y", \
     "plots/MIN3P_Data_Paper/120Y.txt" using 1:3	title "MIN3P 120Y"


reset
set terminal postscript enhanced
set size square 0.65,0.65
set output "Calcite.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Calcite Volume Fraction"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65



plot "postProcessing/graphCell/3.149928e+08/line.xy" using 1:3 with lines lw 5 title "10 Y", \
     "postProcessing/graphCell/3.1499928e+09/line.xy" using 1:3 with lines lw 5 title "100 Y", \
     "postProcessing/graphCell/3.78e+09/line.xy" using 1:3 with lines lw 5 title "120 Y", \
     "plots/MIN3P_Data_Paper/10Y.txt"  using 1:2	title "MIN3P 10Y", \
     "plots/MIN3P_Data_Paper/100Y.txt" using 1:2	title "MIN3P 100Y", \
     "plots/MIN3P_Data_Paper/120Y.txt" using 1:2	title "MIN3P 120Y"


reset
set terminal postscript enhanced
set size square 0.65,0.65
set output "HydraulicHead.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Hydraulic head (m)"
##set yrange [0:0.7]
##set ytics 0.1
set key at 1.0,0.0025 


g = 9.81

plot "postProcessing/graphCell/3.149928e+08/line.xy" using 1:($4/g) with lines lw 5 title "10 Y", \
     "postProcessing/graphCell/3.1499928e+09/line.xy" using 1:($4/g) with lines lw 5 title "100 Y", \
     "postProcessing/graphCell/3.78e+09/line.xy" using 1:($4/g) with lines lw 5 title "120 Y", \
     "plots/MIN3P_Data_Paper/10Y.txt"  using 1:4	title "MIN3P 10Y", \
     "plots/MIN3P_Data_Paper/100Y.txt" using 1:4	title "MIN3P 100Y", \
     "plots/MIN3P_Data_Paper/120Y.txt" using 1:4	title "MIN3P 120Y"



reset
set terminal postscript enhanced
set size square 0.65,0.65
set output "Flux.eps"
##set xrange [0 : 1]
set xlabel "Times (Y)"
set ylabel "Flux (m3/s)"
##set yrange [0:0.7]
##set ytics 0.1
set key at 1.0,0.0025 



plot "postProcessing/patchFlowRate(patch=outlet)/0/surfaceFieldValue.dat" using ($1/(3600*24*365.25)):($2*(3600*24)/0.01) with lines lw 5 title "10 Y", \
     "plots/MIN3P_Data_Paper/Flux.txt"  using 1:2	title "MIN3P"

