Example files are [here](YFeO3)

This example studies YFeO3, calculating S(Q,E), convoluting with a resolution function, and integrating over axes directions. 
This is published in [Phys. Rev. B 89, 014420](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.014420)

CommonFunctions.cpp
```cpp
// required imports
#include <cmath>
#include <string>
#include "CommonFunctions.h"

// avoid appending SpinWaveGenie:: to everything in the namespace
using namespace SpinWaveGenie;

SpinWaveGenie::SpinWave createModel()
{
    // Cell object stores the basis vectors and Sublattice objects
    Cell cell;
    // set a,b,c,alpha,beta,gamma in units of Angstroms,Degrees
    cell.setBasisVectors(5.4,5.4,7.63675323681,90.0,90.0,90.0);
    
    // parameters used for calculating angles
    double delta = 0.00516;
    double phi = 0.0032;
    
    Sublattice Fe1;
    //each Sublattice contains:
    //  a unique name,
    std::string name1 = "Fe1";
    Fe1.setName(name1);
    //  type (used to calculate the magnetic form factor)
    Fe1.setType("FE3");
    //  magnetic moment (magnitude,theta (radians), phi (radians))
    Fe1.setMoment(2.5,M_PI/2.0-delta,M_PI+phi);
    //add sublattice to Cell cell.
    cell.addSublattice(Fe1);
    // add atom to sublattice name1 at position (0,0.5,0) in reduced lattice units.
    cell.addAtom(name1,0.0,0.5,0.0);

    //repeat for Sublattice Fe2,Fe3 and Fe4
    Sublattice Fe2;
    std::string name2 = "Fe2";
    Fe2.setName(name2);
    Fe2.setType("FE3");
    Fe2.setMoment(2.5,M_PI/2.0-delta,phi);
    cell.addSublattice(Fe2);
    cell.addAtom(name2,0.0,0.5,0.5);

    Sublattice Fe3;
    std::string name3 = "Fe3";
    Fe3.setName(name3);
    Fe3.setType("FE3");
    Fe3.setMoment(2.5,M_PI/2.0-delta,M_PI-phi);
    cell.addSublattice(Fe3);
    cell.addAtom(name3,0.5,0.0,0.5);

    Sublattice Fe4;
    std::string name4 = "Fe4";
    Fe4.setName(name4);
    Fe4.setType("FE3");
    Fe4.setMoment(2.5,M_PI/2.0-delta,2.0*M_PI-phi);
    cell.addSublattice(Fe4);
    cell.addAtom(name4,0.5,0.0,0.0);
    
    // add Cell cell to the builder
    SpinWaveBuilder builder(cell);
    
    // defining interactions.
    InteractionFactory interactions;
   
    Vector3 xhat(1.0,0.0,0.0); 
    // get anisotropy interaction named "Ka" with value "-0.0055" in the "xhat" direction on sublattice name1.
    // pass this interaction to the builder.
    // repeat for interactions name2,name3,name4 
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.0055,xhat,name1));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.0055,xhat,name2));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.0055,xhat,name3));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.0055,xhat,name4));

    Vector3 zhat(0.0,0.0,1.0);
    // get anisotropy interaction named "Kc" with value "-0.00305" in the "zhat" direction on sublattice name1.
    // pass this interaction to the builder.
    // repeat for interactions name2,name3,name4 
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00305,zhat,name1));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00305,zhat,name2));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00305,zhat,name3));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00305,zhat,name4));

    // get exchange interaction named "J1" with value "-4.77" between sublattices name1 and name2.
    // limit atoms to those between 3.8 and 4.3 angstroms.
    // pass this interaction to the builder.
    // repeat for interactions between name1 and name4, name2 and name3, name3 and name4. 
    builder.addInteraction(interactions.getExchange("J1",-4.77,name1,name2,3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.77,name1,name4,3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.77,name3,name2,3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.77,name3,name4,3.8,4.3));

    // get exchange interaction named "J2" with value "-0.21" between sublattices name1 and name1.
    // limit atoms to those between 5.3 and 5.5 angstroms.
    // pass this interaction to the builder.
    // repeat for interactions between other pairs. 
    builder.addInteraction(interactions.getExchange("J2",-0.21,name1,name1,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name2,name2,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name3,name3,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name4,name4,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name1,name3,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name2,name4,5.3,5.5));

    Vector3 yhat(0.0,1.0,0.0);
    // get Dzyaloshinskii-Moriya interaction named "D1" with value "-0.074" in the yhat direction between sublattices name4 and name1.
    // limit atoms to those between 3.8 and 4.3 angstroms.
    // pass this interaction to the builder.
    // repeat for interactions between other pairs. 
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D1",-0.074,yhat,name4,name1,3.8,4.3));
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D1",-0.074,yhat,name2,name3,3.8,4.3));

    // get Dzyaloshinskii-Moriya interaction named "D2" with value "-0.028" in the zhat direction between sublattices name4 and name1.
    // limit atoms to those between 3.8 and 4.3 angstroms.
    // pass this interaction to the builder.
    // repeat for interactions between other pairs. 
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D2",-0.028,zhat,name4,name1,3.8,4.3));
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D2",-0.028,zhat,name3,name2,3.8,4.3));

    //return SpinWave object containing the cell, sublattices and interactions.
    return builder.createElement();
}
```

TwoDimensionalCut.cpp
```cpp
#include "CommonFunctions.h"

using namespace std;
using namespace SpinWaveGenie;

int main()
{
    SpinWave SW = createModel();
    
    PointsAlongLine Line;
    Line.setFirstPoint(2.0,-1.5,-3.0);
    Line.setFinalPoint(2.0, 1.5,-3.0);
    Line.setNumberPoints(201);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    Energies energies(0.0, 80.0, 201);
    
    TwoDimGaussian resinfo;
    resinfo.a = 1109.0;
    resinfo.b = 0.0;
    resinfo.c = 0.48;
    resinfo.tol = 5.0e-3;
    resinfo.direction = Vector3(0.0,1.0,0.0);
    
    unique_ptr<SpinWavePlot> res(new TwoDimensionResolutionFunction(resinfo, SW, energies));
    
    HKLDirections direction;
    direction.addDirection(0.0,0.0,1.0,0.2);
    direction.addDirection(1.0,0.0,0.0,0.2);
    //direction.addDirection(0.0,1.0,0.0,0.05);
    
    unique_ptr<SpinWavePlot> asdf(new IntegrateAxes(std::move(res),direction,1.0e-2));

    TwoDimensionalCut twodimcut;
    twodimcut.setFilename("YFeO3_2_K_m3");
    twodimcut.setPlotObject(move(asdf));
    twodimcut.setPoints(kPoints);
    twodimcut.save();

    Line.setFirstPoint(2.0,-1.5,-2.0);
    Line.setFinalPoint(2.0, 1.5,-2.0);
    Line.setNumberPoints(201);
    kPoints = Line.getPoints();

    twodimcut.setFilename("YFeO3_2_K_m2");
    twodimcut.setPoints(kPoints);
    twodimcut.save();
    return 0;
}
```
