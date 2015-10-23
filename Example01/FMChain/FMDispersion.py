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

dispersion = np.transpose(np.loadtxt('FMChain.txt'))

ax.plot(dispersion[0],dispersion[3],c='black')
ax.set_ylabel(r'Energy Transfer (meV)')
ax.set_xlabel(r'$(\xi,0,0)$ (rlu)')
ax.set_xlim(left=0.0,right=3.0)
ax.set_ylim(bottom=0.0,top=5.0)
plt.savefig('FMDispersion.png',dpi=400,bbox_inches='tight')
plt.close()
