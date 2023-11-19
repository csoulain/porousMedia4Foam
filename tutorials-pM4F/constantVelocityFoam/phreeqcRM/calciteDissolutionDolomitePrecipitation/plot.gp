
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

plot "postProcessing/graphCell/1200/line.xy" using 1:2 with lines lw 5 title "1200 s", \
     "postProcessing/graphCell/2400/line.xy" using 1:2 with lines lw 5 title "2400 s", \
     "postProcessing/graphCell/3600/line.xy" using 1:2 with lines lw 5 title "3600 s", \
     "plots/CalciteDissDoloPrecip1200.op"  using 1:($10*Vm_calc*poro*1e-3)	title "1200s PHREEQC", \
     "plots/CalciteDissDoloPrecip2400.op" using  1:($10*Vm_calc*poro*1e-3)	title "2400s PHREEQC", \
     "plots/CalciteDissDoloPrecip3600.op" using  1:($10*Vm_calc*poro*1e-3)	title "3600s PHREEQC"






set terminal postscript enhanced
set size square 0.65,0.65
set output "Dolomite.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Dolomite Volume Fraction"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65

Vm_dolo = 64.5
poro = 0.4

plot "postProcessing/graphCell/1200/line.xy" using 1:7 with lines lw 5 title "1200 s", \
     "postProcessing/graphCell/2400/line.xy" using 1:7 with lines lw 5 title "2400 s", \
     "postProcessing/graphCell/3600/line.xy" using 1:7 with lines lw 5 title "3600 s", \
     "plots/CalciteDissDoloPrecip1200.op"  using 1:($8*poro*Vm_dolo*1e-3)	title "1200s PHREEQC", \
     "plots/CalciteDissDoloPrecip2400.op" using 1:($8*poro*Vm_dolo*1e-3)	title "2400s PHREEQC", \
     "plots/CalciteDissDoloPrecip3600.op" using 1:($8*poro*Vm_dolo*1e-3)	title "3600s PHREEQC"



set terminal postscript enhanced
set size square 0.65,0.65
set output "Cl.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Cl concentration"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65

plot "postProcessing/graphCell/1200/line.xy" using 1:3 with lines lw 5 title "1200 s", \
     "postProcessing/graphCell/2400/line.xy" using 1:3 with lines lw 5 title "2400 s", \
     "postProcessing/graphCell/3600/line.xy" using 1:3 with lines lw 5 title "3600 s", \
     "plots/CalciteDissDoloPrecip1200.op"  using 1:4	title "1200s PHREEQC", \
     "plots/CalciteDissDoloPrecip2400.op" using 1:4	title "2400s PHREEQC", \
     "plots/CalciteDissDoloPrecip3600.op" using 1:4	title "3600s PHREEQC"




set terminal postscript enhanced
set size square 0.65,0.65
set output "Mg.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Mg concentration"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65

plot "postProcessing/graphCell/1200/line.xy" using 1:6 with lines lw 5 title "1200 s", \
     "postProcessing/graphCell/2400/line.xy" using 1:6 with lines lw 5 title "2400 s", \
     "postProcessing/graphCell/3600/line.xy" using 1:6 with lines lw 5 title "3600 s", \
     "plots/CalciteDissDoloPrecip1200.op"  using 1:7	title "1200s PHREEQC", \
     "plots/CalciteDissDoloPrecip2400.op" using 1:7	title "2400s PHREEQC", \
     "plots/CalciteDissDoloPrecip3600.op" using 1:7	title "3600s PHREEQC"




set terminal postscript enhanced
set size square 0.65,0.65
set output "pH.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "pH"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65

plot "postProcessing/graphCell/1200/line.xy" using 1:8 with lines lw 5 title "1200 s", \
     "postProcessing/graphCell/2400/line.xy" using 1:8 with lines lw 5 title "2400 s", \
     "postProcessing/graphCell/3600/line.xy" using 1:8 with lines lw 5 title "3600 s", \
     "plots/CalciteDissDoloPrecip1200.op"  using 1:3	title "1200s PHREEQC", \
     "plots/CalciteDissDoloPrecip2400.op" using 1:3	title "2400s PHREEQC", \
     "plots/CalciteDissDoloPrecip3600.op" using 1:3	title "3600s PHREEQC"



set terminal postscript enhanced
set size square 0.65,0.65
set output "Ca.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Ca concentration"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65

plot "postProcessing/graphCell/1200/line.xy" using 1:5 with lines lw 5 title "1200 s", \
     "postProcessing/graphCell/2400/line.xy" using 1:5 with lines lw 5 title "2400 s", \
     "postProcessing/graphCell/3600/line.xy" using 1:5 with lines lw 5 title "3600 s", \
     "plots/CalciteDissDoloPrecip1200.op"  using 1:5	title "1200s PHREEQC", \
     "plots/CalciteDissDoloPrecip2400.op" using 1:5	title "2400s PHREEQC", \
     "plots/CalciteDissDoloPrecip3600.op" using 1:5	title "3600s PHREEQC"



set terminal postscript enhanced
set size square 0.65,0.65
set output "C.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "C concentration"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65

plot "postProcessing/graphCell/1200/line.xy" using 1:4 with lines lw 5 title "1200 s", \
     "postProcessing/graphCell/2400/line.xy" using 1:4 with lines lw 5 title "2400 s", \
     "postProcessing/graphCell/3600/line.xy" using 1:4 with lines lw 5 title "3600 s", \
     "plots/CalciteDissDoloPrecip1200.op"  using 1:6	title "1200s PHREEQC", \
     "plots/CalciteDissDoloPrecip2400.op" using 1:6	title "2400s PHREEQC", \
     "plots/CalciteDissDoloPrecip3600.op" using 1:6	title "3600s PHREEQC"


