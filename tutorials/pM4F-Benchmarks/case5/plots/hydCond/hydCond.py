import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plots/porosity/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(.6,.9,5e-9,1e-7), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Hydraulic conductivity [m/s]')
plt.title('')

#MIN3P - Paper Data
DistMIN = []
HC10MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_10Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       HC10MIN.append(float(row[2]))
plt.plot(DistMIN, HC10MIN,'r-*',label='MIN3P',linewidth=2.5)

HC100MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_100Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       HC100MIN.append(float(row[2]))
plt.plot(DistMIN, HC100MIN,'r-x',label='MIN3P',linewidth=2.5)

HH300MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_300Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       HH300MIN.append(float(row[2]))
plt.plot(DistMIN, HH300MIN,'r-v',label='MIN3P')

#ToughRe
Dist10TR = []
HC10TR = []
with open('plots/Tough_WebExtraction/hyCond/hyCond_10y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       HC10TR.append(float(row[1]))
plt.plot(Dist10TR, HC10TR,'m-*',label='TR',linewidth=2.5)

Dist100TR = []
HC100TR = []
with open('plots/Tough_WebExtraction/hyCond/hyCond_100y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))
       HC100TR.append(float(row[1]))
plt.plot(Dist100TR, HC100TR,'m-x',label='TR',linewidth=2.5)

Dist300TR = []
HC300TR = []
with open('plots/Tough_WebExtraction/hyCond/hyCond_300y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist300TR.append(float(row[0]))
       HC300TR.append(float(row[1]))
plt.plot(Dist300TR, HC300TR,'m-v',label='ToughRe',linewidth=2.5)

#pM4F hyd conductivity results
Dist = []
HC10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       HC10.append(float(row[3])*(1000*9.81)/(0.001))
plt.plot(Dist, HC10,'c-*',label='pM4F',linewidth=2.5)

HC100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       HC100.append(float(row[3])*(1000*9.81)/(0.001))
plt.plot(Dist, HC100,'c-x',label='pM4F',linewidth=2.5)

HC300 = []
with open('postProcessing/singleGraph/9.45e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       HC300.append(float(row[3])*(1000*9.81)/(0.001))
plt.plot(Dist, HC300,'c-v',label='pM4F',linewidth=2.5)


#axes.set_xlim([0.,1])
#axes.set_ylim([1e-12,1e-2])

plt.ylim(top=1e-2)
plt.ylim(bottom=1e-12)
plt.xlim(left=0.)
plt.xlim(right=1.)

plt.yscale('log')

x_ticks=np.arange(0.,1.001,0.2)
#y_ticks=np.arange(0.001,1.001)
plt.xticks(x_ticks)
#plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)

plt.savefig('B3_hydCond.eps', format='eps')
plt.show()
