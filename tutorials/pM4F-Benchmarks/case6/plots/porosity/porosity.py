import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plots/porosity/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(.6,.95,0.002,0.01), zorder=1)

plt.xlabel('Distance [m]', fontsize = 18)
plt.ylabel('Porosity [-]', fontsize = 18)
plt.title('')

#MIN3P - Paper Data
DistMIN = []
Poro10MIN = []
with open('plots/MIN3P_Data_Paper/B5_MIN3P_10Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       Poro10MIN.append(float(row[1]))
plt.plot(DistMIN, Poro10MIN,'r-*',label='MIN3P', linewidth=2.5)

Poro100MIN = []
with open('plots/MIN3P_Data_Paper/B5_MIN3P_100Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Poro100MIN.append(float(row[1]))
plt.plot(DistMIN, Poro100MIN,'r-x',label='MIN3P',linewidth=2.5)

Poro300MIN = []
with open('plots/MIN3P_Data_Paper/B5_MIN3P_300Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Poro300MIN.append(float(row[1]))
plt.plot(DistMIN, Poro300MIN,'r-v',label='MIN3P',linewidth=2.5)

#ToughR
Dist10TR = []
Poro10TR = []
with open('plots/Tough_WebExtraction/Porosity/B5_10y_poro.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       Poro10TR.append(float(row[1]))
plt.plot(Dist10TR, Poro10TR,'m-*',label='ToughReact',linewidth=2.5)

Dist100TR = []
Poro100TR = []
with open('plots/Tough_WebExtraction/Porosity/B5_100y_poro.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))	
       Poro100TR.append(float(row[1]))
plt.plot(Dist100TR, Poro100TR,'m-x',label='Toughreact',linewidth=2.5)

Dist300TR = []
Poro300TR = []
with open('plots/Tough_WebExtraction/Porosity/B5_300y_poro.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist300TR.append(float(row[0]))  
       Poro300TR.append(float(row[1]))
plt.plot(Dist300TR, Poro300TR,'m-v',label='Toughreact',linewidth=2.5)

#pM4F Porosity results
Dist = []
Poro10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Poro10.append(float(row[1]))
plt.plot(Dist, Poro10,'c-*',label='pM4F',linewidth=2.5)

Poro100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro100.append(float(row[1]))
plt.plot(Dist, Poro100,'c-x',label='pM4F',linewidth=2.5)

Poro300 = []
with open('postProcessing/singleGraph/9.45e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro300.append(float(row[1]))
plt.plot(Dist, Poro300,'c-v',label='pM4F',linewidth=2.5)


#axes.set_xlim([0.,1])
#axes.set_ylim([0.00001,1])

plt.ylim(top=1)
plt.ylim(bottom=0.00001)
plt.xlim(left=0.)
plt.xlim(right=1.)

plt.yscale("log")

x_ticks=np.arange(0.,1.001,0.2)
#y_ticks=np.arange(0.001,1.001)
plt.xticks(x_ticks)
#plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)


plt.savefig('B5_porosity.eps', format='eps')
plt.show()
