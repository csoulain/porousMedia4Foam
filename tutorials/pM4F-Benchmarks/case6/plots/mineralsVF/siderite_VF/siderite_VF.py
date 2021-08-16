import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plots/porosity/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(.6,.95,0.1,0.14), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Volume fraction of siderite [-]')
plt.title('')

#MIN3P - Paper Data
DistMIN = []
Sid10MIN = []
with open('plots/MIN3P_Data_Paper/B5_MIN3P_10Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       Sid10MIN.append(float(row[10]))
plt.plot(DistMIN, Sid10MIN,'r-*',label='MIN3P',linewidth=2.5)

Sid100MIN = []
with open('plots/MIN3P_Data_Paper/B5_MIN3P_100Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Sid100MIN.append(float(row[10]))
plt.plot(DistMIN, Sid100MIN,'r-x',label='MIN3P',linewidth=2.5)

Sid300MIN = []
with open('plots/MIN3P_Data_Paper/B5_MIN3P_300Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Sid300MIN.append(float(row[9]))
plt.plot(DistMIN, Sid300MIN,'r-v',label='MIN3P',linewidth=2.5)

#ToughRe
Dist10TR = []
Sid10TR = []
with open('plots/Tough_WebExtraction/minVF/sid_VF/sid_VF_10y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       Sid10TR.append(float(row[1]))
plt.plot(Dist10TR, Sid10TR,'m-*',label='TR',linewidth=2.5)

Dist100TR = []
Sid100TR = []
with open('plots/Tough_WebExtraction/minVF/sid_VF/sid_VF_100y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))
       Sid100TR.append(float(row[1]))
plt.plot(Dist100TR, Sid100TR,'m-x',label='TR',linewidth=2.5)

Dist300TR = []
Sid300TR = []
with open('plots/Tough_WebExtraction/minVF/sid_VF/sid_VF_300y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist300TR.append(float(row[0]))
       Sid300TR.append(float(row[1]))
plt.plot(Dist300TR, Sid300TR,'m-v',label='TR',linewidth=2.5)


#pM4F Siderite VF results
Dist = []
Sid10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Sid10.append(float(row[8]))
plt.plot(Dist, Sid10,'c-*',label='pM4F',linewidth=2.5)

Sid100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Sid100.append(float(row[8]))
plt.plot(Dist, Sid100,'c-x',label='pM4F',linewidth=2.5)

Sid300 = []
with open('postProcessing/singleGraph/9.45e+09/line_eps_p_K_Ys.Calcite_Ys.Gypsum_Ys.Ferrihydrite_Ys.Gibbsite_Ys.Siderite_Ys.Jarosite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Sid300.append(float(row[8]))
plt.plot(Dist, Sid300,'c-v',label='pM4F',linewidth=2.5)

#axes.set_xlim([0.,1])
#axes.set_ylim([0.,0.25])

plt.ylim(top=0.25)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=1.)

x_ticks=np.arange(0.,1.001,0.2)
y_ticks=np.arange(0.,.25001,0.05)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)


plt.savefig('B5_sidVF.eps', format='eps')
plt.show()
