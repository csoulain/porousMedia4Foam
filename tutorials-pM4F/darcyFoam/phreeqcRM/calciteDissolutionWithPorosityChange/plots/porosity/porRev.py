import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

#im = image.imread('plots/porosity/legend.eps')
#fig, ax = plt.subplots()
#ax.imshow(im,aspect='auto', extent=(0.275,0.9,0.37,0.45), zorder=1)

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
#plt.plot(DistMIN, Poro10MIN,'r-*',label='MIN3P',linewidth=5)
##plt.plot(DistMIN, Poro10MIN,'r',marker="*",label='MIN3P',linewidth=3,markersize=9)
plt.plot(DistMIN, Poro10MIN,'r',label='MIN3P',linewidth=3)

Poro100MIN = []
with open('plots/MIN3P_Data_Paper/100Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Poro100MIN.append(float(row[2]))
#plt.plot(DistMIN, Poro100MIN,'r-x',label='MIN3P',linewidth=5)
##plt.plot(DistMIN, Poro100MIN, 'r', marker="x",linewidth=3,markersize=9)
plt.plot(DistMIN, Poro100MIN,'r',linewidth=3)

Poro120MIN = []
with open('plots/MIN3P_Data_Paper/120Y.txt','r') as csvfile:
   plots = csv.reader(csvfile, delimiter='\t')
   for row in plots:
       Poro120MIN.append(float(row[2]))
#plt.plot(DistMIN, Poro120MIN,'r-v',label='MIN3P',linewidth=5)
##plt.plot(DistMIN, Poro120MIN,'r', marker="v",linewidth=3,markersize=9)
plt.plot(DistMIN,Poro120MIN,'r',linewidth=3)

#pM4F Porosity results
Dist = []
Poro10 = []
with open('postProcessing/singleGraph/3.149928e+08/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Poro10.append(float(row[1]))
#plt.plot(Dist, Poro10,'c-*',label='pM4F',linewidth=5)
##plt.plot(Dist, Poro10,'c',marker="*",linewidth=3,markersize=9)
plt.plot(Dist, Poro10,'c',marker="s",markevery=5,markersize=10,linestyle='None')

Poro100 = []
with open('postProcessing/singleGraph/3.1499928e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro100.append(float(row[1]))
#plt.plot(Dist, Poro100,'c-x',label='pM4F',linewidth=5)
##plt.plot(Dist, Poro100,'c',marker="x",linewidth=3,markersize=9)
plt.plot(Dist, Poro100,'c',marker="s",markevery=5,markersize=10,linestyle='None')

Poro120 = []
with open('postProcessing/singleGraph/3.78e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro120.append(float(row[1]))
#plt.plot(Dist, Poro120,'c-v',label='pM4F',linewidth=5)
##plt.plot(Dist, Poro120,'c',marker="v",linewidth=3,markersize=9)
plt.plot(Dist, Poro120,'c',marker="s",markevery=5,markersize=10,linestyle='None')

#pM4F-dbs Porosity results
DistDBS = []
Poro10DBS = []
with open('../B1_dbs/postProcessing/singleGraph/3.149928e+08/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistDBS.append(float(row[0]))
       Poro10DBS.append(float(row[1]))
#plt.plot(DistDBS, Poro10DBS,'k-*',label='pM4F',linewidth=5)
##plt.plot(DistDBS, Poro10DBS,'k',marker="*",linewidth=3,markersize=9)
plt.plot(DistDBS, Poro10DBS,'k',marker="o",markevery=5,markersize=6,linestyle='None')

Poro100DBS = []
with open('../B1_dbs/postProcessing/singleGraph/3.1499928e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro100DBS.append(float(row[1]))
##plt.plot(DistDBS, Poro100DBS,'k',marker="x",linewidth=3,markersize=9)
plt.plot(DistDBS, Poro100DBS,'k',marker="o",markevery=5,markersize=6,linestyle='None')

Poro120DBS = []
with open('../B1_dbs/postProcessing/singleGraph/3.78e+09/line_eps_Ys.Calcite_p.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Poro120DBS.append(float(row[1]))
##plt.plot(DistDBS, Poro120DBS,'k',marker="v",linewidth=3,markersize=9)
plt.plot(DistDBS, Poro120DBS,'k',marker="o",markevery=5,markersize=6,linestyle='None')


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

