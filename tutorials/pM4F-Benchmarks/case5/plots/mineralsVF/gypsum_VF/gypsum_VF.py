import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()
im = image.imread('plots/mineralsVF/gypsum_VF/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(.6,.95,0.35,0.425), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Volume fraction of gypsum [-]')
plt.title('')

#MIN3P - Paper Data
DistMIN = []
Gyp10MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_10Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       Gyp10MIN.append(float(row[6]))
plt.plot(DistMIN, Gyp10MIN,'r-*',label='MIN3P',linewidth=2.5)

Gyp100MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_100Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Gyp100MIN.append(float(row[6]))
plt.plot(DistMIN, Gyp100MIN,'r-x',label='MIN3P',linewidth=2.5)

Gyp300MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_300Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Gyp300MIN.append(float(row[5]))
plt.plot(DistMIN, Gyp300MIN,'r-v',label='MIN3P',linewidth=2.5)

#ToughRe
Dist10TR = []
Gyp10TR = []
with open('plots/Tough_WebExtraction/minVF/gyp_VF/gyp_VF_10y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       Gyp10TR.append(float(row[1]))
plt.plot(Dist10TR, Gyp10TR,'m-*',label='TR',linewidth=2.5)

Dist100TR = []
Gyp100TR = []
with open('plots/Tough_WebExtraction/minVF/gyp_VF/gyp_VF_100y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))
       Gyp100TR.append(float(row[1]))
plt.plot(Dist100TR, Gyp100TR,'m-x',label='TR',linewidth=2.5)

Dist300TR = []
Gyp300TR = []
with open('plots/Tough_WebExtraction/minVF/gyp_VF/gyp_VF_300y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist300TR.append(float(row[0]))
       Gyp300TR.append(float(row[1]))
plt.plot(Dist300TR, Gyp300TR,'m-v',label='TR',linewidth=2.5)

#pM4F Gypsum VF results
Dist = []
Gyp10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Gyp10.append(float(row[5]))
plt.plot(Dist, Gyp10,'c-*',label='pM4F',linewidth=2.5)

Gyp100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Gyp100.append(float(row[5]))
plt.plot(Dist, Gyp100,'c-x',label='pM4F',linewidth=2.5)

Gyp300 = []
with open('postProcessing/singleGraph/9.45e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Gyp300.append(float(row[5]))
plt.plot(Dist, Gyp300,'c-v',label='pM4F', linewidth=2.5)

#axes.set_xlim([0.,1])
#axes.set_ylim([0.,0.5])

plt.ylim(top=0.5)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=1.)

x_ticks=np.arange(0.,1.001,0.2)
y_ticks=np.arange(0.,.5001,0.05)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)


plt.savefig('B3_gypVF.eps', format='eps')
plt.show()
