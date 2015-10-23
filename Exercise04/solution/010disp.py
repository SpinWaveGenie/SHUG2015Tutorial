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

dispersion = np.transpose(np.loadtxt('010disp.txt'))

ax.plot(dispersion[1],dispersion[3],c='black')
ax.plot(dispersion[1],dispersion[4],c='black')
ax.plot(dispersion[1],dispersion[5],c='black')
ax.plot(dispersion[1],dispersion[6],c='black')
ax.plot(dispersion[1],dispersion[7],c='black')
ax.plot(dispersion[1],dispersion[8],c='black')
ax.set_ylabel(r'Energy (meV)')
ax.set_xlabel(r'$(2,\xi,0)$ (rlu)')
ax.set_xlim(left=-6.0,right=4.0)
ax.set_ylim(bottom=0.0,top=30.0)
plt.savefig('along010.png',dpi=400,bbox_inches='tight')
plt.close()
