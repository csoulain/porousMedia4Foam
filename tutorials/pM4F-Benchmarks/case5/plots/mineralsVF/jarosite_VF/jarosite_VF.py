import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plots/mineralsVF/jarosite_VF//legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(.06,.09,0.045,0.055), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Volume fraction of jarosite [-]')
plt.title('')

#MIN3P - Paper Data
DistMIN = []
Jar10MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_10Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       Jar10MIN.append(float(row[8]))
plt.plot(DistMIN, Jar10MIN,'r-*',label='MIN3P',linewidth=2.5)

Jar100MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_100Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Jar100MIN.append(float(row[8]))
plt.plot(DistMIN, Jar100MIN,'r-x',label='MIN3P',linewidth=2.5)

Jar300MIN = []
with open('plots/MIN3P_Data_Paper/B3_MIN3P_300Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Jar300MIN.append(float(row[7]))
plt.plot(DistMIN, Jar300MIN,'r-v',label='MIN3P',linewidth=2.5)

#ToughRe
Dist10TR = []
Jar10TR = []
with open('plots/Tough_WebExtraction/minVF/jar_VF/jar_VF_10y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       Jar10TR.append(float(row[1]))
plt.plot(Dist10TR, Jar10TR,'m-*',label='TR',linewidth=2.5)

Dist100TR = []
Jar100TR = []
with open('plots/Tough_WebExtraction/minVF/jar_VF/jar_VF_100y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))
       Jar100TR.append(float(row[1]))
plt.plot(Dist100TR, Jar100TR,'m-x',label='TR',linewidth=2.5)

Dist300TR = []
Jar300TR = []
with open('plots/Tough_WebExtraction/minVF/jar_VF/jar_VF_300y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist300TR.append(float(row[0]))
       Jar300TR.append(float(row[1]))
plt.plot(Dist300TR, Jar300TR,'m-v',label='TR',linewidth=2.5)

#pM4F Jarosite VF results
Dist = []
Jar10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Jar10.append(float(row[9]))
plt.plot(Dist, Jar10,'c-*',label='pM4F',linewidth=2.5)

Jar100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Jar100.append(float(row[9]))
plt.plot(Dist, Jar100,'c-x',label='pM4F',linewidth=2.5)

Jar300 = []
with open('postProcessing/singleGraph/9.45e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Jar300.append(float(row[9]))
plt.plot(Dist, Jar300,'c-v',label='pM4F',linewidth=2.5)

#axes.set_xlim([0.,.1])
#axes.set_ylim([0.,0.06])

plt.ylim(top=0.06)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=0.1)

x_ticks=np.arange(0.,.1001,0.02)
y_ticks=np.arange(0.,.06001,0.01)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)


plt.savefig('B3_jarVF.eps', format='eps')
plt.show()
