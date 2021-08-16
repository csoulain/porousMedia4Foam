import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plots/calcite_VF/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(1.2,1.8,0.05,0.125), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Volume fraction of calcite [-]')
plt.title('')

#MIN3P - Paper Data
DistMIN = []
Calc10MIN = []
with open('plots/MIN3P_Data_Paper/B2_10y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       Calc10MIN.append(float(row[1]))
plt.plot(DistMIN, Calc10MIN,'r-*',label='MIN3P',linewidth=2.5)

Calc100MIN = []
with open('plots/MIN3P_Data_Paper/B2_100y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Calc100MIN.append(float(row[1]))
plt.plot(DistMIN, Calc100MIN,'r-x',label='MIN3P',linewidth=2.5)

#Calc120MIN = []
#with open('../MIN3P_Data_Paper/120Y.txt','r') as csvfile:
#   plots = csv.reader(csvfile, delimiter='\t')
#   for row in plots:
#       Calc120MIN.append(float(row[1]))
#plt.plot(DistMIN, Calc120MIN,'r-o',label='MIN3P')

#ToughRe
Dist10TR = []
Calc10TR = []
with open('plots/Tough_WebExtraction/calc_VF/calc_VF_10y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       Calc10TR.append(float(row[1]))
plt.plot(Dist10TR, Calc10TR,'m-',label='TR',linewidth=2.5)

Dist100TR = []
Calc100TR = []
with open('plots/Tough_WebExtraction/calc_VF/calc_VF_100y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))
       Calc100TR.append(float(row[1]))
plt.plot(Dist100TR, Calc100TR,'m-',label='TR',linewidth=2.5)

#Calc120MIN = []
#with open('../MIN3P_Data_Paper/120Y.txt','r') as csvfile:
#   plots = csv.reader(csvfile, delimiter='\t')
#   for row in plots:
#       Calc120MIN.append(float(row[1]))
#plt.plot(DistMIN, Calc120MIN,'r-o',label='MIN3P')


#pM4F Calcite VF results
Dist = []
Calc10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_p_Ys.Gypsum_K_Ys.Calcite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Calc10.append(float(row[5]))
plt.plot(Dist, Calc10,'c-*',label='pM4F',linewidth=2.5)

Calc100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_p_Ys.Gypsum_K_Ys.Calcite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Calc100.append(float(row[5]))
plt.plot(Dist, Calc100,'c-x',label='pM4F',linewidth=2.5)

#Calc120 = []
#with open('../../postProcessing/singleGraph/3.78e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
#   plots = csv.reader(csvfile, delimiter=' ')
#   for row in plots:
#       Calc120.append(float(row[2]))
#plt.plot(Dist, Calc120,'k-o',label='pM4F')

#axes.set_xlim([0.,2])
#axes.set_ylim([0.,0.4])

plt.ylim(top=0.4)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=2.)

x_ticks=np.arange(0.,2.001,0.2)
y_ticks=np.arange(0.,.4001,0.1)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)

plt.savefig('B2_calcVF.eps', format='eps')
plt.show()
