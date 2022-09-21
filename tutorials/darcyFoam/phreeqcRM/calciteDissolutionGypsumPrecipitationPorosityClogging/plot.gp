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

##3.149928e+08,3.1499928e+09

plot "plots/MIN3P_Data_Paper/B2_10y.txt"  using 1:4 every 3 pt 1			title "10Y MIN3P", \
     "plots/MIN3P_Data_Paper/B2_100y.txt" using 1:4 every 3 pt 2			title "100Y MIN3P", \
     "plots/Tough_WebExtraction/Porosity/B2_10y_poro.txt" using 1:2	every 3 pt 3 	title "10Y TOUCHREACT", \
     "plots/Tough_WebExtraction/Porosity/B2_100y_poro.txt" using 1:2	every 3 pt 4 	title "100Y TOUCHREACT" ,\
     "postProcessing/graphCell/3.149928e+08/line.xy" using 1:2 with lines  title "10 Y", \
     "postProcessing/graphCell/3.1499928e+09/line.xy" using 1:2 with lines title "100 Y"




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


Vm_calc = 36.9
poro = 0.4


plot "plots/MIN3P_Data_Paper/B2_10y.txt"  using 1:2 every 3 pt 1			title "10Y MIN3P", \
     "plots/MIN3P_Data_Paper/B2_100y.txt" using 1:2 every 3 pt 2			title "100Y MIN3P", \
     "plots/Tough_WebExtraction/calc_VF/calc_VF_10y.txt" using 1:2	every 3 pt 3 	title "10Y TOUCHREACT", \
     "plots/Tough_WebExtraction/calc_VF/calc_VF_100y.txt" using 1:2	every 3 pt 4 	title "100Y TOUCHREACT" ,\
     "postProcessing/graphCell/3.149928e+08/line.xy" using 1:6 with lines  title "10 Y", \
     "postProcessing/graphCell/3.1499928e+09/line.xy" using 1:6 with lines title "100 Y"




set terminal postscript enhanced
set size square 0.65,0.65
set output "Gypsum.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Gypsum Volume Fraction"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65


Vm_gyps = 74.2
Vm_dolo = 64.5
poro = 0.4

plot "plots/MIN3P_Data_Paper/B2_10y.txt"  using 1:3 every 3 pt 1			title "10Y MIN3P", \
     "plots/MIN3P_Data_Paper/B2_100y.txt" using 1:3 every 3 pt 2			title "100Y MIN3P", \
     "plots/Tough_WebExtraction/gyp_VF/gyp_VF_10y.txt" using 1:2	every 3 pt 3 	title "10Y TOUCHREACT", \
     "plots/Tough_WebExtraction/gyp_VF/gyp_VF_100y.txt" using 1:2	every 3 pt 4 	title "100Y TOUCHREACT" ,\
     "postProcessing/graphCell/3.149928e+08/line.xy" using 1:4 with lines  title "10 Y", \
     "postProcessing/graphCell/3.1499928e+09/line.xy" using 1:4 with lines title "100 Y"



set terminal postscript enhanced
set size square 0.65,0.65
set output "HydraulicHead.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Hydraulic head (m)"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65

rho = 1e3
g   = 9.81

plot "plots/MIN3P_Data_Paper/B2_10y.txt"  using 1:6 every 3 pt 1			title "10Y MIN3P", \
     "plots/MIN3P_Data_Paper/B2_100y.txt" using 1:6 every 3 pt 2			title "100Y MIN3P", \
     "plots/Tough_WebExtraction/hyHead/hh_10Y.txt" using 1:2	every 3 pt 3 	title "10Y TOUCHREACT", \
     "plots/Tough_WebExtraction/hyHead/hh_100Y.txt" using 1:2	every 3 pt 4 	title "100Y TOUCHREACT" ,\
     "postProcessing/graphCell/3.149928e+08/line.xy"  using 1:($3/rho/g) with lines  title "10 Y", \
     "postProcessing/graphCell/3.1499928e+09/line.xy" using 1:($3/rho/g) with lines title "100 Y"




set terminal postscript enhanced
set size square 0.65,0.65
set output "HydraulicConductivity.eps"
##set xrange [0 : 1]
set logscale y
set xlabel "Distance (m)"
set ylabel "Hydraulic conductivity (m/s)"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65


rho = 1e3
g   = 9.81
mu = 1e-3


plot "plots/MIN3P_Data_Paper/B2_10y.txt"  using 1:5 every 3 pt 1			title "10Y MIN3P", \
     "plots/MIN3P_Data_Paper/B2_100y.txt" using 1:5 every 3 pt 2			title "100Y MIN3P", \
     "plots/Tough_WebExtraction/hyCond/hyCond_10y.txt" using 1:2	every 3 pt 3 	title "10Y TOUCHREACT", \
     "plots/Tough_WebExtraction/hyCond/hyCond_100y.txt" using 1:2	every 3 pt 4 	title "100Y TOUCHREACT" ,\
     "postProcessing/graphCell/3.149928e+08/line.xy"  using 1:($5*rho*g/mu) with lines  title "10 Y", \
     "postProcessing/graphCell/3.1499928e+09/line.xy" using 1:($5*rho*g/mu) with lines title "100 Y"



reset
set terminal postscript enhanced
set size square 0.65,0.65
set output "Flux.eps"
set xrange [0 : 150]
set xlabel "Times (Y)"
set ylabel "Flux (m3/s)"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 1.0,0.0025 


plot "postProcessing/patchFlowRate(patch=outlet)/0/surfaceFieldValue.dat" using ($1/(3600*24*365.25)):($2*(3600*24)/0.01) with lines lw 5 title "OpenFOAM-PHREEQC", \
     "plots/MIN3P_Data_Paper/B2_Flux_150y.txt"  using 1:2 every 15 pt 1	title "MIN3P" ,\
     "plots/Tough_WebExtraction/Flux.txt"       using 1:2 pt 4 title "TOUGHREACT"


