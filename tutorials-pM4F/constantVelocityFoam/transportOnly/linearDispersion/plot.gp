reset
##set terminal postscript enhanced
set terminal postscript eps enhanced color font 'Verdana,16'
set size square 0.65,0.65
set output "breakthroughCurve.eps"
##set xrange [0 : 1]
set xlabel "Time (s)"
set ylabel "Concentration"
##set yrange [0:0.7]
##set ytics 0.1
##set key at 0.80,0.65

eps = 0.5 
L   = 10.   ## m
U   = 0.2    ## m/s
D   = 1.e-5 ## m^2/s
C0  = 1.    ## concentration

U1 = U/eps

alphaL = 0.5 ## m

n = 1

Deff=(eps**n)*(D+alphaL*U)


C(x,t)=0.5*C0*(erfc((L-U1*t)*0.5/(sqrt(Deff*t))) +  exp((U*L)/(Deff))*erfc((L+U1*t)*0.5/(sqrt(Deff*t))))


##plot [t=0:60] C(L,t) with lines lw 5 title "Analytical"




##plot [t=0:60] A(L,t) with lines lw 5 title "A(x,t)", \
##	      B(L,t) with lines lw 5 title "B(x,t)"


##plot [x=-100:100] erfc(x) with lines lw 5 title "Analytical"


plot [t=0:100] C(L,t) with lines lc 1 lw 5 title "Analytical", \
	"postProcessing/probes(points=((1000)),fields=(Ci))/0/Ci" using 1:2 every 3 with points  pt 82 lc 1 title "Numerical"

