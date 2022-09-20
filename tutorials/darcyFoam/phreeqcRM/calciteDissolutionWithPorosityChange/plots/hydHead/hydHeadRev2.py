import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

#im = image.imread('plots/hydHead/legend.eps')
#fig, ax = plt.subplots()
#ax.imshow(im,aspect='auto', extent=(0.2,0.8,0.001,0.0025), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Hydraulic head [m]')
plt.title('')

#ToughRe
Dist10TR = []
HH10TR = []
with open('plots/Tough_Data_WebExtract/hyHead/10Y_hh.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       HH10TR.append(float(row[1]))
plt.plot(Dist10TR, HH10TR,'c-',linewidth=2.5)

Dist100TR = []
HH100TR = []
with open('plots/Tough_Data_WebExtract/hyHead/100Y_hh.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))
       HH100TR.append(float(row[1]))
plt.plot(Dist100TR, HH100TR,'c-',linewidth=2.5)

Dist120TR = []
HH120TR = []
with open('plots/Tough_Data_WebExtract/hyHead/120Y_hh.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist120TR.append(float(row[0]))  
       HH120TR.append(float(row[1]))
plt.plot(Dist120TR, HH120TR,'c:',linewidth=2.5)

#MIN3P - Paper Data
DistMIN = []
HH10MIN = []
with open('plots/MIN3P_Data_Paper/10Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       HH10MIN.append(float(row[3]))
plt.plot(DistMIN, HH10MIN,'r-',linewidth=2.5,marker="o",markevery=5,markersize=6)

HH100MIN = []
with open('plots/MIN3P_Data_Paper/100Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       HH100MIN.append(float(row[3]))
plt.plot(DistMIN, HH100MIN,'r-.',linewidth=2.5,marker="o",markevery=5,markersize=6)

HH120MIN = []
with open('plots/MIN3P_Data_Paper/120Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       HH120MIN.append(float(row[3]))
plt.plot(DistMIN, HH120MIN,'r:',linewidth=2.5,marker="o",markevery=5,markersize=6)

#pM4F Porosity results
Dist = []
HH10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       HH10.append(float(row[3])/(1000*9.81))
plt.plot(Dist, HH10,'k-',marker="s",markevery=5,markersize=6)

HH100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       HH100.append(float(row[3])/(1000*9.81))
plt.plot(Dist, HH100,'k-.',marker="s",markevery=5,markersize=6)

HH120 = []
with open('postProcessing/singleGraph/3.78e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       HH120.append(float(row[3])/(1000*9.81))
plt.plot(Dist, HH120,'k:',marker="s",markevery=5,markersize=6)

#axes.set_xlim([0.,2])
#axes.set_ylim([0.,0.008])

plt.ylim(top=0.008)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=2.)


x_ticks=np.arange(0.,2.001,0.5)
y_ticks=np.arange(0.,0.008001,0.002)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.minorticks_on()
#plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
#plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='out',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='out',length=7.5, width=3)

plt.savefig('B1_hydHead.eps', format='eps')
plt.show()
