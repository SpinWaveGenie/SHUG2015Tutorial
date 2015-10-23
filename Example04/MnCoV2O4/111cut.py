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

x_axis = np.loadtxt('111cut1.x')
E_array = np.loadtxt('111cut1.y')
image1 = np.loadtxt('111cut1.mat')
image2 = np.loadtxt('111cut2.mat')
image3 = np.loadtxt('111cut3.mat')
image_array = image1+image2+image3

im = ax.pcolorfast(x_axis[:,2],E_array,image_array,vmin=0.0,vmax=0.4)
ax.set_ylabel(r'Energy (meV)')
ax.set_xlabel(r'distance $(\xi,\xi,\xi)$')
ax.set_xlim(left=0.5,right=1.5)
ax.set_ylim(bottom=0.0,top=15.0)
fig.colorbar(im);
plt.savefig('111cut.png',dpi=400,bbox_inches='tight')
plt.close()
