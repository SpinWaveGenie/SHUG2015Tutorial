#include <cmath>
#include <string>
#include "SpinWaveGenie/SpinWaveGenie.h"
#include "nlopt.hpp"

using namespace std;
using namespace SpinWaveGenie;

SpinWaveBuilder getBuilder(double theta_a, double theta_b, double gamma)
{
    double eta,J;
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

double myfunc(const std::vector<double> &x, std::vector<double> & /*grad*/, void * data)
{
    double *gamma = reinterpret_cast<double*>(data);
    SpinWaveBuilder builder = getBuilder(x[0],x[1],*gamma);
    return builder.getEnergy();
}

std::vector<double> findAngles(double gamma)
{
    
    double minf = 0.0;
    std::vector<double> lb = {0.0,0.0};
    std::vector<double> ub = {M_PI,M_PI};
    std::vector<double> thetaphi = {M_PI/4.0,M_PI/4.0};
    
    nlopt::opt opt(nlopt::LN_COBYLA,2);
    opt.set_upper_bounds(ub);
    opt.set_lower_bounds(lb);
    
    opt.set_ftol_abs(1.0e-13);
    opt.set_maxeval(50000000);
    
    opt.set_min_objective(myfunc,&gamma);
    opt.optimize(thetaphi,minf);
    
    return thetaphi;
}

int main()
{
    
    
    for(int i = 0; i<5;++i)
    {
      vector<double> x = findAngles(static_cast<double>(i));
      cout << x[0] << " " << x[1]-M_PI_2 << endl;
    
      SpinWaveBuilder builder = getBuilder(x[0], x[1],static_cast<double>(i));
      SpinWave test = builder.createElement();
    
      PointsAlongLine Line;
      Line.setFirstPoint(1.0,0.0,0.0);
      Line.setFinalPoint(0.0,0.0,0.0);
      Line.setNumberPoints(101);
      ThreeVectors<double> points100 = Line.getPoints();
    
      Line.setFirstPoint(0.0,0.0,0.0);
      Line.setFinalPoint(0.0,1.0,0.0);
      Line.setNumberPoints(101);
      ThreeVectors<double> points010 = Line.getPoints();
    
      SpinWaveDispersion dispersion;
      dispersion.setFilename("SC1Chain_"+std::to_string(i)+".txt");
      dispersion.setGenie(test);
      dispersion.setPoints(points100);
      dispersion.setPoints(points010);
    
      dispersion.save();
    }
    return 0;
}
