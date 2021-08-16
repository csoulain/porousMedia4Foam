import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('calciteGrainData.txt')

x = data[:,0]
y = data[:,1]
plt.plot(x,y,'r-',linewidth=3)

z = data[:,2]
plt.plot(x,z,'k-',linewidth=3)

axes = plt.gca()
axes.set_xlim([0,2700])
axes.set_ylim([0,3.5e-8])

plt.show()
