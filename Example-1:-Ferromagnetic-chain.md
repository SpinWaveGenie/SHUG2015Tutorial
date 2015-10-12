example files are at 
https://github.com/SpinWaveGenie/SpinWaveGenie/tree/master/examples/FMChain

Dispersion.cpp
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "SpinWaveGenie/SpinWaveGenie.h"

using namespace std;
using namespace SpinWaveGenie;

int main()
{
    Cell cell;
    cell.setBasisVectors(1.0,10.0,10.0,90.0,90.0,90.0);
    
    Sublattice spin0;
    string name0 = "Spin0";
    spin0.setName(name0);
    spin0.setType("NONE");
    spin0.setMoment(1.0,0.0,0.0);
    cell.addSublattice(spin0);
    cell.addAtom(name0,0.0,0.0,0.0);

    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    builder.addInteraction(interactions.getExchange("J",1.0,name0,name0,0.9,1.1));

    SpinWave test = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(0.0,0.0,0.0);
    Line.setFinalPoint(3.0,0.0,0.0);
    Line.setNumberPoints(61);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    SpinWaveDispersion dispersion;
    dispersion.setFilename("FMChain.txt");
    dispersion.setGenie(test);
    dispersion.setPoints(kPoints);
    
    dispersion.save();
    return 0;
}
```

TwoDimensionalCut.cpp
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "SpinWaveGenie/SpinWaveGenie.h"
using namespace std;
using namespace SpinWaveGenie;

int main()
{
    Cell cell;
    cell.setBasisVectors(1.0,10.0,10.0,90.0,90.0,90.0);
    
    Sublattice spin0;
    string name0 = "Spin0";
    spin0.setName(name0);
    spin0.setType("NONE");
    spin0.setMoment(1.0,0.0,0.0);
    cell.addSublattice(spin0);
    cell.addAtom(name0,0.0,0.0,0.0);

    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    builder.addInteraction(interactions.getExchange("J",1.0,name0,name0,0.9,1.1));

    SpinWave SW = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(0.0,0.0,0.0);
    Line.setFinalPoint(3.0,0.0,0.0);
    Line.setNumberPoints(401);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    Energies energies(0.0, 5.0, 401);
    
    OneDimensionalFactory factory;
    auto gauss = factory.getGaussian(0.25,1.0e-5);
    
    unique_ptr<SpinWavePlot> res(new EnergyResolutionFunction(move(gauss), SW,energies));
    
    TwoDimensionalCut twodimcut;
    twodimcut.setFilename("FMcut");
    twodimcut.setPlotObject(move(res));
    twodimcut.setPoints(kPoints);
    twodimcut.save();
    return 0;
}
```
FM.py
```python
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
image_array = np.loadtxt('FMcut.mat')
E_array = np.loadtxt('FMcut.y')
x_axis = np.loadtxt('FMcut.x')

print len(image_array)
print len(image_array[0])
print np.max(image_array)
im = ax.pcolorfast(x_axis[:,0],E_array,image_array,vmin=0.0,vmax=1.0)
ax.scatter(dispersion[0],dispersion[3],c='white')
ax.set_ylabel(r'Energy (meV)')
ax.set_xlabel(r'distance $(\xi,0,0)$')
ax.set_xlim(left=0.0,right=3.0)
ax.set_ylim(bottom=0.0,top=5.0)
fig.colorbar(im);
plt.savefig('FMChain.png',dpi=400,bbox_inches='tight')
plt.close()
```
<img src="https://github.com/SpinWaveGenie/SpinWaveGenie/blob/master/examples/FMChain/FMChain.png" width="auto" height="600px" />

PowderAverage.cpp
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "SpinWaveGenie/SpinWaveGenie.h"

using namespace std;
using namespace SpinWaveGenie;

int main()
{
    Cell cell;
    cell.setBasisVectors(1.0,10.0,10.0,90.0,90.0,90.0);
    Sublattice spin0;
    string name0 = "Spin0";
    spin0.setName(name0);
    spin0.setType("NONE");
    spin0.setMoment(1.0,0.0,0.0);
    cell.addSublattice(spin0);
    cell.addAtom(name0,0.0,0.0,0.0);

    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    builder.addInteraction(interactions.getExchange("J",1.0,name0,name0,0.9,1.1));

    SpinWave SW = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(0.0,0.0,0.0);
    Line.setFinalPoint(0.0,0.0,3.0*2.0*M_PI);
    Line.setNumberPoints(401);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    Energies energies(0.0, 5.0, 401);
    
    OneDimensionalFactory factory;
    auto gauss = factory.getGaussian(0.25,1.0e-2);
    
    unique_ptr<SpinWavePlot> res(new EnergyResolutionFunction(move(gauss), SW,energies));
    unique_ptr<SpinWavePlot> cut(new IntegrateThetaPhi(move(res),2.0e-2));

    TwoDimensionalCut twodimcut;
    twodimcut.setFilename("FMPowderAverage");
    twodimcut.setPlotObject(move(cut));
    twodimcut.setPoints(kPoints);
    twodimcut.save();
    return 0;
}
```
FMavg.py
```python
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

dispersion = np.loadtxt('FMChain.txt')
x_axis = np.loadtxt('FMPowderAverage.x')
E_array = np.loadtxt('FMPowderAverage.y')
image_array = np.loadtxt('FMPowderAverage.mat')

print len(image_array)
print len(image_array[0])
print np.max(image_array)
im = ax.pcolorfast(x_axis[:,2]/(2.0*np.pi),E_array,image_array,vmin=0.0,vmax=0.3)
ax.scatter(dispersion[:,0],dispersion[:,3],c='white')
ax.set_ylabel(r'Energy (meV)')
ax.set_xlabel(r'distance $(\xi,0,0)$')
ax.set_xlim(left=0.0,right=3.0)
ax.set_ylim(bottom=0.0,top=5.0)
fig.colorbar(im);
plt.savefig('FMavg.png',dpi=400,bbox_inches='tight')
plt.close()
```
<img src="https://github.com/SpinWaveGenie/SpinWaveGenie/blob/master/examples/FMChain/FMavg.png" width="auto" height="600px" />