reset
set terminal postscript enhanced
set size square 0.65,0.65
set output "CalciteVolumeFractionThermEq.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Calcite volume fraction"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65



plot "postProcessing/graphCell/1200/line.xy" using 1:2 with lines lw 5 title "20 min", \
     "postProcessing/graphCell/2400/line.xy" using 1:2 with lines lw 5 title "40 min", \
     "postProcessing/graphCell/3600/line.xy" using 1:2 with lines lw 5 title "60 min", \
     "phreeqcDataExtraction/T20" using 1:7	title "Phreeqc", \
     "phreeqcDataExtraction/T40" using 1:7	title "Phreeqc", \
     "phreeqcDataExtraction/T60" using 1:7	title "Phreeqc"
     



reset
set terminal postscript enhanced
set size square 0.65,0.65
set output "ClIonConcentrationThermEq.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Cl- Concentration"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65



plot "postProcessing/graphCell/1200/line.xy" using 1:3 with lines lw 5 title "20 min Cl-", \
     "phreeqcDataExtraction/T20" using 1:4	title "Phreeqc 20 min Cl-", \
     "postProcessing/graphCell/2400/line.xy" using 1:3 with lines lw 5 title "40 min Cl-", \
     "phreeqcDataExtraction/T40" using 1:4	title "Phreeqc 40 min Cl-", \
     "postProcessing/graphCell/3600/line.xy" using 1:3 with lines lw 5 title "60 min Cl-", \
     "phreeqcDataExtraction/T60" using 1:4	title "Phreeqc 60 min Cl-" 


reset
set terminal postscript enhanced
set size square 0.65,0.65
set output "CaIonConcentrationThermEq.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Ca2+ Concentration"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65



plot "postProcessing/graphCell/1200/line.xy" using 1:5 with lines lw 5 title "20 min Ca2+", \
     "phreeqcDataExtraction/T20" using 1:5	title "Phreeqc 20 min Ca2+", \
     "postProcessing/graphCell/2400/line.xy" using 1:5 with lines lw 5 title "40 min Ca2+", \
     "phreeqcDataExtraction/T40" using 1:5	title "Phreeqc 40 min Ca2+", \
     "postProcessing/graphCell/3600/line.xy" using 1:5 with lines lw 5 title "60 min Ca2+", \
     "phreeqcDataExtraction/T60" using 1:5	title "Phreeqc 60 min Ca2+" 
