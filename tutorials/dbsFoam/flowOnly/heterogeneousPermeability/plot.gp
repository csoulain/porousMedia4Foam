reset 

k1 = 1e-11
k2 = 5e-11
k3 = 1e-12
k4 = 5e-11

l = 0.5
L = 4.

S = 0.1*0.1

K = L/(l/k1+l/k2+l/k3+l/k4)
mu = 1.e-3
rho = 1.e3


P0 = 1.e5/rho
P4 = 0


V(x) = -K/mu*(P4-P0)/L*rho
print "V = ", V(0)

P3 = P4 + K/L*l/k4*(P0-P4)
P2 = P3 + k4/l*l/k3*(P3-P4)
P1 = P2 + k3/l*l/k2*(P2-P3)

print P3
print P2
print P1


f1(x) = x < l ? P0-(P0-P1)/l*x :  1/0
f2(x) = x < l ? 1/0 : x<2*l ? (2*P1-P2)-(P1-P2)/l*x : 1/0
f3(x) = x < 2*l ? 1/0 : x<3*l ? (P2+(P2-P3)*2)-(P2-P3)/l*x : 1/0
f4(x) = x < 3*l ? 1/0 : x<4*l ? (P3+(P3-P4)*3)-(P3-P4)/l*x  : 1/0

f(x) = x < l ? P0-(P0-P1)/l*x : x<2*l ? (2*P1-P2)-(P1-P2)/l*x : x<3*l ? (P2+(P2-P3)*2)-(P2-P3)/l*x : x<4*l ? (P3+(P3-P4)*3)-(P3-P4)/l*x  : 1/0



set terminal postscript enhanced
set size square 0.65,0.65
set output "Permeability.eps"
##set yrange [1e-12 : 5e-11]
set xlabel "Distance (m)"
set ylabel "Permeability"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65



plot "postProcessing/graphCell/0/line.xy" using 1:6 with lines lw 5 title "darcyFoam"


set terminal postscript enhanced
set size square 0.65,0.65
set output "Velocity.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Velocity (m/s)"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65



plot "postProcessing/graphCell/100/line.xy" using 1:3 with lines lw 5 title "darcyFoam",\
	V(x) with  linespoints pointinterval 3 ps 2 title "analytical"




set terminal postscript enhanced
set size square 0.65,0.65
set output "Pressure.eps"
##set xrange [0 : 1]
set xlabel "Distance (m)"
set ylabel "Pressure (Pa)"
##set yrange [0:0.7]
##set ytics 0.1



plot "postProcessing/graphCell/100/line.xy" using 1:2 with lines lw 5 title "darcyFoam",\
     f(x)  with  linespoints pointinterval 3 ps 2 title "analytical"



