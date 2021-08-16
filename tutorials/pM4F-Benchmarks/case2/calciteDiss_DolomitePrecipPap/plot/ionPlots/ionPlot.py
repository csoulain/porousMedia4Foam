import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plot/ionPlots/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(0.275,0.475,0.0175,0.0275), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Ion concentration [mol/l]')
plt.title('')

Dist = []
Ca20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Ca20min.append(float(row[4]))
plt.plot(Dist, Ca20min,'k*',label='pM4F')

Dist = []
Ca40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Ca40min.append(float(row[4]))
plt.plot(Dist, Ca40min,'kx',label='pM4F')

Dist = []
Ca60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
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
       DistPh.append(float(row[0]))
       Ca60minPh.append(float(row[4]))
plt.plot(DistPh, Ca60minPh,'ro',label='Phreeqc')

###Carbonate
Dist = []
C20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C20min.append(float(row[3]))
plt.plot(Dist, C20min,'g*',label='pM4F')

Dist = []
C40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C40min.append(float(row[3]))
plt.plot(Dist, C40min,'gx',label='pM4F')

Dist = []
C60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C60min.append(float(row[3]))
plt.plot(Dist, C60min,'go',label='pM4F')


###Phreeqc-Carbonate
DistPh = []
C20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       C20minPh.append(float(row[5]))
plt.plot(DistPh, C20minPh,'b*',label='Phreeqc')

C40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       C40minPh.append(float(row[5]))
plt.plot(DistPh, C40minPh,'bx',label='Phreeqc')

C60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       C60minPh.append(float(row[5]))
plt.plot(DistPh, C60minPh,'bo',label='Phreeqc')

###Chloride
Dist = []
Cl20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Cl20min.append(float(row[2]))
plt.plot(Dist, Cl20min,'k*',label='pM4F')

Dist = []
Cl40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Cl40min.append(float(row[2]))
plt.plot(Dist, Cl40min,'kx',label='pM4F')

Dist = []
Cl60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
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
plt.plot(DistPh, Cl20minPh,'r*',label='Phreeqc')

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

###pM4F-Mg
Dist = []
Mg20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Mg20min.append(float(row[5]))
plt.plot(Dist, Mg20min,'C1*',label='pM4F')

Mg40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Mg40min.append(float(row[5]))
plt.plot(Dist, Mg40min,'C1x',label='pM4F')

Mg60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Mg60min.append(float(row[5]))
plt.plot(Dist, Mg60min,'C1o',label='pM4F')

###Phreeqc-Mg
DistPh = []
Mg20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       Mg20minPh.append(float(row[6]))
plt.plot(DistPh, Mg20minPh,'m*',label='Phreeqc')

Mg40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Mg40minPh.append(float(row[6]))
plt.plot(DistPh, Mg40minPh,'mx',label='Phreeqc')

Mg60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Mg60minPh.append(float(row[6]))
plt.plot(DistPh, Mg60minPh,'mo',label='Phreeqc')

plt.ylim(top=0.035)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=0.5)

ax.annotate("$Cl^{-}$",
            xy=(0.027, 0.031), xycoords='data',
            xytext=(0.02, 0.026), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

plt.savefig('ionConPlot_CalcDolo.eps', format='eps')
plt.show()

#Copy of above data for generating a zoomed plot for concentrations
im = image.imread('plot/ionPlots/legend.eps')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(0.0125,0.145,0.011,0.0145), zorder=1)

plt.xlabel('Distance [m]')
plt.ylabel('Ion concentration [mol/l]')
plt.title('')

Dist = []
Ca20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Ca20min.append(float(row[4]))
plt.plot(Dist, Ca20min,'k*',label='pM4F')

Dist = []
Ca40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Ca40min.append(float(row[4]))
plt.plot(Dist, Ca40min,'kx',label='pM4F')

Dist = []
Ca60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
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
       DistPh.append(float(row[0]))
       Ca60minPh.append(float(row[4]))
plt.plot(DistPh, Ca60minPh,'ro',label='Phreeqc')

###Carbonate
Dist = []
C20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C20min.append(float(row[3]))
plt.plot(Dist, C20min,'g*',label='pM4F')

Dist = []
C40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C40min.append(float(row[3]))
plt.plot(Dist, C40min,'gx',label='pM4F')

Dist = []
C60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       C60min.append(float(row[3]))
plt.plot(Dist, C60min,'go',label='pM4F')


###Phreeqc-Carbonate
DistPh = []
C20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       C20minPh.append(float(row[5]))
plt.plot(DistPh, C20minPh,'b*',label='Phreeqc')

C40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       C40minPh.append(float(row[5]))
plt.plot(DistPh, C40minPh,'bx',label='Phreeqc')

C60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       C60minPh.append(float(row[5]))
plt.plot(DistPh, C60minPh,'bo',label='Phreeqc')

###Chloride
Dist = []
Cl20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Cl20min.append(float(row[2]))
plt.plot(Dist, Cl20min,'k*',label='pM4F')

Dist = []
Cl40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Cl40min.append(float(row[2]))
plt.plot(Dist, Cl40min,'kx',label='pM4F')

Dist = []
Cl60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
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
plt.plot(DistPh, Cl20minPh,'r*',label='Phreeqc')

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

###pM4F-Mg
Dist = []
Mg20min = []

with open('dt0_75/postProcessing/singleGraph/1200/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Dist.append(float(row[0]))
       Mg20min.append(float(row[5]))
plt.plot(Dist, Mg20min,'C1*',label='pM4F')

Mg40min = []

with open('dt0_75/postProcessing/singleGraph/2400/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Mg40min.append(float(row[5]))
plt.plot(Dist, Mg40min,'C1x',label='pM4F')

Mg60min = []

with open('dt0_75/postProcessing/singleGraph/3600/line_Ys.Calcite_Y.Cl_Y.C_Y.Ca_Y.Mg_Ys.Dolomite.xy','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Mg60min.append(float(row[5]))
plt.plot(Dist, Mg60min,'C1o',label='pM4F')

###Phreeqc-Mg
DistPh = []
Mg20minPh = []

with open('plot/phreeqcDataExtraction/T20','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       DistPh.append(float(row[0]))
       Mg20minPh.append(float(row[6]))
plt.plot(DistPh, Mg20minPh,'m*',label='Phreeqc')

Mg40minPh = []

with open('plot/phreeqcDataExtraction/T40','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Mg40minPh.append(float(row[6]))
plt.plot(DistPh, Mg40minPh,'mx',label='Phreeqc')

Mg60minPh = []

with open('plot/phreeqcDataExtraction/T60','r') as csvfile:
   plots = csv.reader(csvfile, delimiter=' ')
   for row in plots:
       Mg60minPh.append(float(row[6]))
plt.plot(DistPh, Mg60minPh,'mo',label='Phreeqc')

plt.ylim(top=0.015)
plt.ylim(bottom=0.)
plt.xlim(left=0.)
plt.xlim(right=0.5)

ax.annotate("$Ca^{+2}$",
            xy=(0.45, 0.006), xycoords='data',
            xytext=(0.45, 0.004), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

ax.annotate("$Mg^{+2}$",
            xy=(0.38, 0.013), xycoords='data',
            xytext=(0.3, 0.013), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

ax.annotate("$CO_{3}^{-2}$",
            xy=(0.15, 0.0085), xycoords='data',
            xytext=(0.08, 0.0085), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

plt.savefig('ionConPlot_CalcDoloZoom.eps', format='eps')
plt.show()
