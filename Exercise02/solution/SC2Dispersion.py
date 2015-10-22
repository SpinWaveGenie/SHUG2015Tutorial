import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt

rc('font',**{'family':'serif','serif':['Times'],'size':16.0})
rc('text', usetex=True)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(5.0*2,3.0*2)
ax = plt.subplot(1,2,1)

disp_1 = np.transpose(np.loadtxt('SC2Chain_0.txt'))
disp_2 = np.transpose(np.loadtxt('SC2Chain_4.txt'))
disp_3 = np.transpose(np.loadtxt('SC2Chain_8.txt'))
disp_4 = np.transpose(np.loadtxt('SC2Chain_10.txt'))

ax.plot(disp_1[0,0:100],disp_1[3,0:100],c='black')
ax.plot(disp_1[0,0:100],disp_1[4,0:100],c='black')
ax.plot(disp_2[0,0:100],disp_2[3,0:100],c='red',ls=':')
ax.plot(disp_2[0,0:100],disp_2[4,0:100],c='red',ls=':')
ax.plot(disp_3[0,0:100],disp_3[3,0:100],c='blue',ls='--')
ax.plot(disp_3[0,0:100],disp_3[4,0:100],c='blue',ls='--')
ax.plot(disp_4[0,0:100],disp_4[3,0:100],c='green',ls='-.')
ax.plot(disp_4[0,0:100],disp_4[4,0:100],c='green',ls='-.')

ax.set_ylabel(r'Energy (meV)')
ax.set_xlabel(r'$(\xi,0,0)$ (rlu)')
ax.set_xlim(left=0.5,right=0.0)
ax.set_ylim(bottom=0.0,top=20.0)

ax2 = plt.subplot(1,2,2)

ax2.plot(disp_1[1,100:],disp_1[3,100:],c='black',label=r'B = 0')
ax2.plot(disp_1[1,100:],disp_1[4,100:],c='black')
ax2.plot(disp_2[1,100:],disp_2[3,100:],c='red',ls=':',label=r'B = 4')
ax2.plot(disp_2[1,100:],disp_2[4,100:],c='red',ls=':')
ax2.plot(disp_3[1,100:],disp_3[3,100:],c='blue',ls='--',label=r'B = 8')
ax2.plot(disp_3[1,100:],disp_3[4,100:],c='blue',ls='--')
ax2.plot(disp_4[1,100:],disp_4[3,100:],c='green',ls='-.',label=r'B = 10')
ax2.plot(disp_4[1,100:],disp_4[4,100:],c='green',ls='-.')
ax2.set_ylabel(r'Energy (meV)')
ax2.set_xlabel(r'$(0,\xi,0)$ (rlu)')
ax2.set_xlim(left=0.0,right=1.0)
ax2.set_ylim(bottom=0.0,top=20.0)
ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('SC2Dispersion.png',dpi=400,bbox_inches='tight')
plt.close()
