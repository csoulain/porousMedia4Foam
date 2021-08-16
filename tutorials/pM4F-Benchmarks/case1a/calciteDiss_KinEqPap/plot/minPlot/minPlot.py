import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plot/minPlot/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(0.365,0.49,0.001,0.0045), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Calcite volume fraction [-]')
plt.title('')

Dist = []
Calcite20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Calcite20min.append(float(row[1]))
plt.plot(Dist, Calcite20min,'k-*',label='pM4F')

Calcite40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite40min.append(float(row[1]))
plt.plot(Dist, Calcite40min,'k-x',label='pM4F')

Calcite60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite60min.append(float(row[1]))
plt.plot(Dist, Calcite60min,'k-o',label='pM4F')

#Phreeqc data
Calcite20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite20minPh.append(float(row[6]))
plt.plot(Dist, Calcite20minPh,'r-*',label='Phreeqc')

Calcite40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite40minPh.append(float(row[6]))
plt.plot(Dist, Calcite40minPh,'r-x',label='Phreeqc')

Calcite60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite60minPh.append(float(row[6]))
plt.plot(Dist, Calcite60minPh,'r-o',label='Phreeqc')

#axes.set_xlim([0.,0.5])
#axes.set_ylim([0.,0.0305])

plt.ylim(top=0.0305)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=0.5)

plt.savefig('minConPlot_CalcDissKin.eps', format='eps')
plt.show()
