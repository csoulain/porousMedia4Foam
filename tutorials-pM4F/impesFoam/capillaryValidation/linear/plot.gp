set terminal postscript enhanced
set size square 0.65,0.65
set xtics 0.2 nomirror

set output "PC_VanGenuchten.eps"
set ytics 10 nomirror
set xlabel "Saturation"
set ylabel "dS/dx (m^{-1})"

pc0 = 0
pcmax = 19620

rhoa = 1e0
rhob = 1e3
gy = -9.81

plot [0.001:1.02] -(rhob-rhoa)*gy/(pcmax-pc0) ti "Analytical",\
     "postProcessing/graphCell/80000/line.xy" using 2:(-$4) ti "Numerical" , \
     "postProcessing/graphCell/10000/line.xy" using 2:(-$4) ti "Numerical = 10000" , \
     "postProcessing/graphCell/40000/line.xy" using 2:(-$4) ti "Numerical = 40000" 

set output "Sb_VanGenuchten.eps"
set ytics 0.2 nomirror
set yrange [0:1.05]
set xlabel "Position (m)"
set ylabel "Saturation"
plot "postProcessing/graphCell/80000/line.xy" using 1:2 w l lt 3 lw 2 ti ""
