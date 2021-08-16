import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

im = image.imread('plots/hydCond/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(1.2,1.8,5e-10,8e-9), zorder=-1)

#axes = plt.gca()

plt.xlabel('Distance [m]')
plt.ylabel('Hydraulic conductivity [m/s]')
plt.title('')

#MIN3P - Paper Data
DistMIN = []
HC10MIN = []
with open('plots/MIN3P_Data_Paper/B2_10y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       HC10MIN.append(float(row[4]))
plt.plot(DistMIN, HC10MIN,'r-*',label='MIN3P',linewidth=2.5)

HC100MIN = []
with open('plots/MIN3P_Data_Paper/B2_100y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       HC100MIN.append(float(row[4]))
plt.plot(DistMIN, HC100MIN,'r-x',label='MIN3P',linewidth=2.5)

#HH120MIN = []
#with open('../MIN3P_Data_Paper/120Y.txt','r') as csvfile:
#   plots = csv.reader(csvfile, delimiter='\t')
#   for row in plots:
#       HH120MIN.append(float(row[3]))
#plt.plot(DistMIN, HH120MIN,'r-o',label='MIN3P')

#ToughRe
Dist10TR = []
HC10TR = []
with open('plots/Tough_WebExtraction/hyCond/hyCond_10y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       HC10TR.append(float(row[1]))
plt.plot(Dist10TR, HC10TR,'m-',label='TR',linewidth=2.5)

Dist100TR = []
HC100TR = []
with open('plots/Tough_WebExtraction/hyCond/hyCond_100y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))
       HC100TR.append(float(row[1]))
plt.plot(Dist100TR, HC100TR,'m-',label='TR',linewidth=2.5)

#HH120MIN = []
#with open('../MIN3P_Data_Paper/120Y.txt','r') as csvfile:
#   plots = csv.reader(csvfile, delimiter='\t')
#   for row in plots:
#       HH120MIN.append(float(row[3]))
#plt.plot(DistMIN, HH120MIN,'r-o',label='MIN3P')

#pM4F hyd conductivity results
Dist = []
HC10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_p_Ys.Gypsum_K_Ys.Calcite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       HC10.append(float(row[4])*(1000*9.81)/(0.001))
plt.plot(Dist, HC10,'c-*',label='pM4F',linewidth=2.5)

HC100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_p_Ys.Gypsum_K_Ys.Calcite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       HC100.append(float(row[4])*(1000*9.81)/(0.001))
plt.plot(Dist, HC100,'c-x',label='pM4F',linewidth=2.5)

#HH120 = []
#with open('../../postProcessing/singleGraph/3.78e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
#   plots = csv.reader(csvfile, delimiter=' ')
#   for row in plots:
#       HH120.append(float(row[3])/(1000*9.81))
#plt.plot(Dist, HH120,'k-o',label='pM4F')


#axes.set_xlim([0.,2])
#axes.set_ylim([1e-12,1e-2])

plt.ylim(top=1e-2)
plt.ylim(bottom=1e-12)
plt.xlim(left=0.)
plt.xlim(right=2.)

plt.yscale('log')

x_ticks=np.arange(0.,2.001,0.2)
#y_ticks=np.arange(0.001,1.001)
plt.xticks(x_ticks)
#plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)

plt.savefig('B2_hydCond.eps', format='eps')
plt.show()
