import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.image as image

#axes = plt.gca()

im = image.imread('plot/minPlot/legend.png')
fig, ax = plt.subplots()
ax.imshow(im,aspect='auto', extent=(3,5,2,3), zorder=-1)

plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.title('')

x1, y1 = [-1,12],[1,4]
x2, y2 = [1,10],[3,2]
plt.plot(x1,y1,x2,y2,marker='o')

axes.set_xlim([0.,0.5])
axes.set_ylim([0.,0.0305])

plt.show()
