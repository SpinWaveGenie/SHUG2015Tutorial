example files at https://github.com/SpinWaveGenie/SpinWaveGenie/tree/master/examples/Villain

reference: http://iopscience.iop.org/0953-8984/21/21/216001

FMDispersion.cpp
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include "SpinWaveGenie/SpinWaveGenie.h"

using namespace std;
using namespace SpinWaveGenie;

int main()
{
    Cell cell;
    cell.setBasisVectors(1.0,2.0,10.0,90.0,90.0,90.0);
    
    Sublattice a1;
    string name1 = "a1";
    a1.setName(name1);
    a1.setType("NONE");
    a1.setMoment(1.0,0.0,0.0);
    cell.addSublattice(a1);
    cell.addAtom(name1,0.0,0.0,0.0);
    
    Sublattice b1;
    string name3 = "b1";
    b1.setName(name3);
    b1.setType("NONE");
    b1.setMoment(1.0,0.0,0.0);
    cell.addSublattice(b1);
    cell.addAtom(name3,0.0,0.5,0.0);
    
    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    double gamma,eta,J,B;
    
    gamma = 3.0;
    eta = -5.0;
    J = 1.0;
    B = 1.0;
    
    builder.addInteraction(interactions.getExchange("J",J,name1,name1,0.8,1.2));
    builder.addInteraction(interactions.getExchange("metaJ",-eta*J,name3,name3,0.8,1.2));
    builder.addInteraction(interactions.getExchange("gammaJ",gamma*J,name1,name3,0.8,1.2));
    
    //Vector3 zhat(0.0,0.0,1.0);
    //builder.addInteraction(interactions.getMagneticField("B",B/2.0,zhat,name1));
    //builder.addInteraction(interactions.getMagneticField("B",B/2.0,zhat,name3));


    SpinWave test = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(1.0,0.0,0.0);
    Line.setFinalPoint(2.0,0.0,0.0);
    Line.setNumberPoints(11);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    SpinWaveDispersion dispersion;
    dispersion.setFilename("AFMChain.txt");
    dispersion.setGenie(test);
    dispersion.setPoints(kPoints);
    
    dispersion.save();
    
    
    double kx(0.0),ky(0.0);
    double R1K = sqrt(pow(eta+1,2)*pow(1-cos(kx),2)+4.0*pow(gamma*cos(ky),2));
    double term = B + 2*gamma + (eta-1.0)*(cos(kx)-1.0);
    cout << term+R1K << " " << term-R1K << endl;
    cout << (R1K-2*gamma*cos(ky))/R1K << " " << (R1K+2*gamma*cos(ky))/R1K << endl;

    
    kx = M_PI;
    R1K = sqrt(pow(eta+1,2)*pow(1-cos(kx),2)+4.0*pow(gamma*cos(ky),2));
    term = B + 2*gamma + (eta-1.0)*(cos(kx)-1.0);
    cout << term+R1K << " " << term-R1K << endl;
    cout << (R1K-2*gamma*cos(ky))/R1K << " " << (R1K+2*gamma*cos(ky))/R1K << endl;

    
    kx = 2.0*M_PI;
    R1K = sqrt(pow(eta+1,2)*pow(1-cos(kx),2)+4.0*pow(gamma*cos(ky),2));
    term = B + 2*gamma + (eta-1.0)*(cos(kx)-1.0);
    cout << term+R1K << " " << term-R1K << endl;
    cout << (R1K-2*gamma*cos(ky))/R1K << " " << (R1K+2*gamma*cos(ky))/R1K << endl;

    return 0;
}
```

SC1Dispersion.cpp
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "SpinWaveGenie/SpinWaveGenie.h"
#include "nlopt.hpp"

using namespace std;
using namespace SpinWaveGenie;

SpinWaveBuilder getBuilder(double theta_a, double theta_b)
{
    double gamma,eta,J;
    
    gamma = 1.0;
    eta = 4.0;
    J = 1.0;
    
    Cell cell;
    cell.setBasisVectors(2.0,2.0,10.0,90.0,90.0,90.0);
    
    Sublattice a1;
    string name_a1 = "a1";
    a1.setName(name_a1);
    a1.setType("NONE");
    a1.setMoment(1.0,theta_a,0.0);
    cell.addSublattice(a1);
    cell.addAtom(name_a1,0.5,0.0,0.0);
    
    Sublattice a2;
    string name_a2 = "a2";
    a2.setName(name_a2);
    a2.setType("NONE");
    a2.setMoment(1.0,theta_a,M_PI);
    cell.addSublattice(a2);
    cell.addAtom(name_a2,0.0,0.0,0.0);
    
    Sublattice b1;
    string name_b1 = "b1";
    b1.setName(name_b1);
    b1.setType("NONE");
    b1.setMoment(1.0,theta_b,0.0);
    cell.addSublattice(b1);
    cell.addAtom(name_b1,0.5,0.5,0.0);
    
    Sublattice b2;
    string name_b2 = "b2";
    b2.setName(name_b2);
    b2.setType("NONE");
    b2.setMoment(1.0,theta_b,M_PI);
    cell.addSublattice(b2);
    cell.addAtom(name_b2,0.0,0.5,0.0);
    
    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    builder.addInteraction(interactions.getExchange("J",J,name_a1,name_a2,0.9,1.1));
    builder.addInteraction(interactions.getExchange("metaJ",-1.0*eta*J,name_b1,name_b2,0.9,1.1));
    builder.addInteraction(interactions.getExchange("gammaJ",gamma*J,name_a1,name_b1,0.9,1.1));
    builder.addInteraction(interactions.getExchange("gammaJ",gamma*J,name_a2,name_b2,0.9,1.1));
    
    return builder;
    
}

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    SpinWaveBuilder builder = getBuilder(x[0],x[1]);
    return builder.getEnergy();
}

std::vector<double> findAngles()
{
    
    double minf = 0.0;
    std::vector<double> lb = {0.0,0.0};
    std::vector<double> ub = {M_PI,M_PI};
    std::vector<double> thetaphi = {M_PI/4.0,M_PI/4.0};
    
    nlopt::opt opt(nlopt::LN_COBYLA,2);
    opt.set_upper_bounds(ub);
    opt.set_lower_bounds(lb);
    
    opt.set_ftol_abs(1.0e-13);
    opt.set_maxeval(5000);
    
    opt.set_min_objective(myfunc,NULL);
    opt.optimize(thetaphi,minf);
    
    return thetaphi;
}

int main()
{
    
    vector<double> x = findAngles();
    cout << x[0] << " " << x[1] << endl;
    
    SpinWaveBuilder builder = getBuilder(x[0], x[1]);
    SpinWave test = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(0.0,0.0,0.0);
    Line.setFinalPoint(0.0,1.0,0.0);
    Line.setNumberPoints(101);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    SpinWaveDispersion dispersion;
    dispersion.setFilename("SC1Chain.txt");
    dispersion.setGenie(test);
    dispersion.setPoints(kPoints);
    
    dispersion.save();
    

    
    return 0;
}
```

SC2Dispersion.cpp
```cpp
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "SpinWaveGenie/Containers/Containers.h"
#include "SpinWaveGenie/Genie/Genie.h"
#include "SpinWaveGenie/Interactions/Interactions.h"
#include "SpinWaveGenie/Plot/Plot.h"

using namespace std;
using namespace SpinWaveGenie;

int main()
{
    
    double gamma,eta,J,B;
    
    gamma = -2.0;
    eta = -2.0;
    J = 1.0;
    B = 8.0;
    
    double theta = M_PI_2;
    theta = acos(-B/(4.0*gamma));
    cout << "Theta = " << theta << endl;
    
    Cell cell;
    cell.setBasisVectors(1.0,2.0,10.0,90.0,90.0,90.0);
    
    Sublattice a1;
    string name1 = "a1";
    a1.setName(name1);
    a1.setType("NONE");
    a1.setMoment(1.0,theta,0.0);
    cell.addSublattice(a1);
    cell.addAtom(name1,0.0,0.0,0.0);
    
    Sublattice b1;
    string name2 = "b1";
    b1.setName(name2);
    b1.setType("NONE");
    b1.setMoment(1.0,theta,M_PI);
    cell.addSublattice(b1);
    cell.addAtom(name2,0.0,0.5,0.0);

    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    builder.addInteraction(interactions.getExchange("J",J,name1,name1,0.9,1.1));
    builder.addInteraction(interactions.getExchange("metaJ",-1.0*eta*J,name2,name2,0.9,1.1));
    builder.addInteraction(interactions.getExchange("gammaJ",gamma*J,name1,name2,0.9,1.1));
    
    Vector3 zhat(0.0,0.0,1.0);
    builder.addInteraction(interactions.getMagneticField("B",B,zhat,name1));
    builder.addInteraction(interactions.getMagneticField("B",B,zhat,name2));

    SpinWave test = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(0.0,1.0,0.0);
    Line.setFinalPoint(0.0,2.0,0.0);
    Line.setNumberPoints(11);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    SpinWaveDispersion dispersion;
    dispersion.setFilename("AFMChain.txt");
    dispersion.setGenie(test);
    dispersion.setPoints(kPoints);
    
    dispersion.save();
    
    return 0;
}