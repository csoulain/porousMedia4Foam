import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plots/porosity/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(0.275,0.9,0.37,0.45), zorder=1)

plt.xlabel('Distance [m]', fontsize = 18)
plt.ylabel('Porosity [-]', fontsize = 18)
plt.title('')

#MIN3P - Paper Data
DistMIN = []
Poro10MIN = []
with open('plots/MIN3P_Data_Paper/10Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       DistMIN.append(float(row[0]))
       Poro10MIN.append(float(row[2]))
plt.plot(DistMIN, Poro10MIN,'r-*',label='MIN3P',linewidth=2.5)

Poro100MIN = []
with open('plots/MIN3P_Data_Paper/100Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Poro100MIN.append(float(row[2]))
plt.plot(DistMIN, Poro100MIN,'r-x',label='MIN3P',linewidth=2.5)

Poro120MIN = []
with open('plots/MIN3P_Data_Paper/120Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Poro120MIN.append(float(row[2]))
plt.plot(DistMIN, Poro120MIN,'r-v',label='MIN3P',linewidth=2.5)

#ToughRe - Web extraction from paper data
Dist10TR = []
Poro10TR = []
with open('plots/Tough_Data_WebExtract/Porosity/10Y_poro.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist10TR.append(float(row[0]))
       Poro10TR.append(float(row[1]))
plt.plot(Dist10TR, Poro10TR,'m-*',label='Toughreact',linewidth=2.5)

Dist100TR = []
Poro100TR = []
with open('plots/Tough_Data_WebExtract/Porosity/100Y_poro.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist100TR.append(float(row[0]))
       Poro100TR.append(float(row[1]))
plt.plot(Dist100TR, Poro100TR,'m-x',label='Toughreact',linewidth=2.5)

Dist120TR = []
Poro120TR = []
with open('plots/Tough_Data_WebExtract/Porosity/120Y_poro.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Dist120TR.append(float(row[0]))
       Poro120TR.append(float(row[1]))
plt.plot(Dist120TR, Poro120TR,'m-v',label='Toughreact',linewidth=2.5)

#pM4F Porosity results
Dist = []
Poro10 = []
with open('../darcyFoam/postProcessing/singleGraph/3.149928e+08/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Poro10.append(float(row[1]))
plt.plot(Dist, Poro10,'c-*',label='pM4F',linewidth=2.5)

Poro100 = []
with open('../darcyFoam/postProcessing/singleGraph/3.1499928e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro100.append(float(row[1]))
plt.plot(Dist, Poro100,'c-x',label='pM4F',linewidth=2.5)

Poro120 = []
with open('../darcyFoam/postProcessing/singleGraph/3.78e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro120.append(float(row[1]))
plt.plot(Dist, Poro120,'c-v',label='pM4F',linewidth=2.5)

#dbs Porosity results
Dist = []
Poro10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Poro10.append(float(row[1]))
plt.plot(Dist, Poro10,'k-*',label='pM4F',linewidth=2.5)

Poro100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro100.append(float(row[1]))
plt.plot(Dist, Poro100,'k-x',label='pM4F',linewidth=2.5)

Poro120 = []
with open('postProcessing/singleGraph/3.78e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro120.append(float(row[1]))
plt.plot(Dist, Poro120,'k-v',label='pM4F',linewidth=2.5)


plt.ylim(top=0.7)
plt.ylim(bottom=0.3)
plt.xlim(left=0.)
plt.xlim(right=2.)

#axes.set_xlim([0.,2])
#axes.set_ylim([0.3,0.7])
#plt.tick_params(labelsize='large')
x_ticks=np.arange(0.,2.001,0.5)
y_ticks=np.arange(0.3,0.7001,0.1)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.minorticks_on()
plt.tick_params(axis='x',which='minor',direction='in',length=5,width=1.5)
plt.tick_params(axis='y',which='minor',direction='in',length=5, width=1.5)
plt.tick_params(axis='x',which='major',direction='in',length=7.5,width=3)
plt.tick_params(axis='y',which='major',direction='in',length=7.5, width=3)

#plt.legend(bbox_to_anchor=(1.005,1))
plt.savefig('B1_porosity.eps', format='eps')
plt.show()

