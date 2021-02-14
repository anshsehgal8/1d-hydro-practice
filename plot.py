import numpy as np 
import matplotlib.pyplot as plt 

dat = np.loadtxt('solution.dat')
x = dat[:,0]
d = dat[:,1]
v = dat[:,2]

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=[6,6])

ax1.plot(x, d, ls='', marker='o', ms=3, color='b')
ax2.plot(x, v, ls='', marker='o', ms=3, color='g')

ax1.set_ylabel('density' , fontsize=11)
ax2.set_ylabel('velocity', fontsize=11)
ax2.set_xlabel('x', fontsize=11)
plt.show()