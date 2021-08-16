import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plot/minPlot/legend.png')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(0.01,0.19,0.024,0.0295), zorder=-1)

plt.xlabel('Distance [m]')
plt.ylabel('Calcite and Dolomite volume fraction [-]')
plt.title('')

#Calcite - pM4F
Dist = []
Calcite20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Calcite20min.append(float(row[1]))
plt.plot(Dist, Calcite20min,'k-*',label='pM4F-Calcite (T=20 min)')

Calcite40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite40min.append(float(row[1]))
plt.plot(Dist, Calcite40min,'k-x',label='pM4F-Calcite (T=40 min)')

Calcite60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite60min.append(float(row[1]))
plt.plot(Dist, Calcite60min,'k-o',label='pM4F-Calcite (T=60 min)')

#Dolomite
Dolomite20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dolomite20min.append(float(row[6]))
plt.plot(Dist, Dolomite20min,'c-*',label='pM4F-Dolomite(T=20 min)')

Dolomite40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dolomite40min.append(float(row[6]))
plt.plot(Dist, Dolomite40min,'c-x',label='pM4F-Dolomite(T=40 min)')

Dolomite60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dolomite60min.append(float(row[6]))
plt.plot(Dist, Dolomite60min,'c-o',label='pM4F-Dolomite(T=40 min)')


#Phreeqc data
#Calcite
Calcite20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite20minPh.append(float(row[8]))
plt.plot(Dist, Calcite20minPh,'r-*',label='Phreeqc-Calcite(T=20 min)')

Calcite40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite40minPh.append(float(row[8]))
plt.plot(Dist, Calcite40minPh,'r-x',label='Phreeqc-Calcite(T=40 min)')

Calcite60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calcite60minPh.append(float(row[8]))
plt.plot(Dist, Calcite60minPh,'r-o',label='Phreeqc-Calcite(T=60 min)')

#Dolomite
Dolomite20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dolomite20minPh.append(float(row[7]))
plt.plot(Dist, Dolomite20minPh,'m-*',label='Phreeqc-Dolomite(T=20 min)')

Dolomite40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dolomite40minPh.append(float(row[7]))
plt.plot(Dist, Dolomite40minPh,'m-x',label='Phreeqc-Dolomite(T=40 min)')

Dolomite60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dolomite60minPh.append(float(row[7]))
plt.plot(Dist, Dolomite60minPh,'m-o',label='Phreeqc-Dolomiteite(T=60 min)')

plt.ylim(top=0.0305)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=0.5)

#axes.set_xlim([0.,0.5])
#axes.set_ylim([0.,0.0305])

#plt.legend(loc='lower center',bbox_to_anchor=(0.5,-0.25),ncol=4,prop={"size":8})
plt.savefig('minCon_CalcDolo.eps', format='eps')
plt.show()
