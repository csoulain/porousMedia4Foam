import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plot/ionPlots/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(0.37,0.49,0.0005,0.002), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Ion concentration [mol/l]')
plt.title('')

Dist = []
Ca20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Ca20min.append(float(row[4]))
plt.plot(Dist, Ca20min,'k*',label='pM4F')

Dist = []
Ca40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Ca40min.append(float(row[4]))
plt.plot(Dist, Ca40min,'kx',label='pM4F')

Dist = []
Ca60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Ca60min.append(float(row[4]))
plt.plot(Dist, Ca60min,'ko',label='pM4F')

##--Phreeqc-Ca
DistPh = []
Ca20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       Ca20minPh.append(float(row[4]))
plt.plot(DistPh, Ca20minPh,'r*',label='Phreeqc')

DistPh = []
Ca40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       Ca40minPh.append(float(row[4]))
plt.plot(DistPh, Ca40minPh,'rx',label='Phreeqc')

DistPh = []
Ca60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Ca60min.append(float(row[4]))
plt.plot(DistPh, Ca60minPh,'ro',label='Phreeqc')

###Carbonate
Dist = []
C20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C20min.append(float(row[3]))
plt.plot(Dist, C20min,'k*',label='pM4F')

Dist = []
C40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C40min.append(float(row[3]))
plt.plot(Dist, C40min,'kx',label='pM4F')

Dist = []
C60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C60min.append(float(row[3]))
plt.plot(Dist, C60min,'ko',label='pM4F')
###Phreeqc-Carbonate
DistPh = []
C20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       C20minPh.append(float(row[5]))
plt.plot(DistPh, C20minPh,'r*',label='Phreeqc')

DistPh = []
C40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       C40minPh.append(float(row[5]))
plt.plot(DistPh, C40minPh,'rx',label='Phreeqc')

DistPh = []
C60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       C60minPh.append(float(row[5]))
plt.plot(DistPh, C60minPh,'ro',label='Phreeqc')

###Chloride
Dist = []
Cl20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Cl20min.append(float(row[2]))
plt.plot(Dist, C20min,'k*',label='pM4F')

Dist = []
Cl40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Cl40min.append(float(row[2]))
plt.plot(Dist, Cl40min,'kx',label='pM4F')

Dist = []
Cl60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Cl60min.append(float(row[2]))
plt.plot(Dist, Cl60min,'ko',label='pM4F')
###Phreeqc-Chloride
DistPh = []
Cl20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       Cl20minPh.append(float(row[3]))
plt.plot(DistPh, C20minPh,'r*',label='Phreeqc')

DistPh = []
Cl40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       Cl40minPh.append(float(row[3]))
plt.plot(DistPh, Cl40minPh,'rx',label='Phreeqc')

DistPh = []
Cl60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       Cl60minPh.append(float(row[3]))
plt.plot(DistPh, Cl60minPh,'ro',label='Phreeqc')

plt.ylim(top=0.012)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=0.5)


#plt.gca().set_xlim([0.,0.5])
#plt.gca().set_ylim([0.,0.012])

ax.annotate("$Cl^{-}$",
            xy=(0.027, 0.0108), xycoords='data',
            xytext=(0.02, 0.009), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

ax.annotate("$Ca^{+2}, CO_{3}^{-2}$",
            xy=(0.2, 0.0082), xycoords='data',
            xytext=(0.2, 0.0095), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )


plt.savefig('ionConPlot_CalcDissKin.eps', format='eps')
plt.show()
