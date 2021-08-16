import matplotlib.pyplot as plt
import numpy as np
import csv

axes = plt.gca()

plt.xlabel('Distance [m]', fontsize = 18)
plt.ylabel('Porosity [-]', fontsize = 18)
plt.title('')

#MIN3P - Paper Data
TimeMIN = []
P1MIN = []
with open('plots/porosity/Min3p_p1','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       TimeMIN.append(float(row[0]))
       P1MIN.append(float(row[1]))
plt.plot(TimeMIN, P1MIN,'r',label='MIN3P', linewidth=5)

TimeMIN1 = []
P2MIN = []
with open('plots/porosity/Min3p_p2','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       TimeMIN1.append(float(row[0]))
       P2MIN.append(float(row[1]))
plt.plot(TimeMIN1, P2MIN,'r-.',label='MIN3P',linewidth=5)


#ToughR
TimeTR = []
P1TR = []
with open('plots/porosity/TR_p1','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       TimeTR.append(float(row[0]))
       P1TR.append(float(row[1]))
plt.plot(TimeTR, P1TR,'m',label='TR', linewidth=5)

TimeTR1 = []
P2TR = []
with open('plots/porosity/TR_p2','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       TimeTR1.append(float(row[0]))
       P2TR.append(float(row[1]))
plt.plot(TimeTR1, P2TR,'m-.',label='TR',linewidth=5)

#pM4F Porosity results
Time = []
P1 = []
with open('postProcessing/probes/4.98e+06/P1','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Time.append(float(row[0])/(86400*365))
       P1.append(float(row[1]))
plt.plot(Time, P1,'c',label='pM4F',linewidth=5)

P2 = []
with open('postProcessing/probes/4.98e+06/P2','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       P2.append(float(row[1]))
plt.plot(Time, P2,'c-.',label='pM4F',linewidth=5)

axes.set_xlim([0.3,300])
axes.set_ylim([0.0,0.6])

#x_ticks=np.arange(0.3,300.001)
y_ticks=np.arange(0.,.6001,0.1)
plt.xscale("log")

#plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)


plt.savefig('B7_porosity.eps', format='eps')
plt.show()
