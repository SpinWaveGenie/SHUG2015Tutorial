import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt

rc('font',**{'family':'serif','serif':['Times'],'size':16.0})
rc('text', usetex=True)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(4.0*2,3.0*2)
ax = plt.subplot(1,1,1)

dispersion = np.transpose(np.loadtxt('111disp.txt'))

ax.plot(dispersion[2],dispersion[3],c='black')
ax.plot(dispersion[2],dispersion[4],c='black')
ax.plot(dispersion[2],dispersion[5],c='black')
ax.plot(dispersion[2],dispersion[6],c='black')
ax.plot(dispersion[2],dispersion[7],c='black')
ax.plot(dispersion[2],dispersion[8],c='black')
ax.set_ylabel(r'Energy (meV)')
ax.set_xlabel(r'$(2+\xi,2+\xi,\xi)$ (rlu)')
ax.set_xlim(left=-2.0,right=1.0)
ax.set_ylim(bottom=0.0,top=30.0)
plt.savefig('111disp.png',dpi=400,bbox_inches='tight')
plt.close()
