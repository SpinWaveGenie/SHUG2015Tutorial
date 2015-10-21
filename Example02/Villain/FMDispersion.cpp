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
    
    Vector3 zhat(0.0,0.0,1.0);
    builder.addInteraction(interactions.getMagneticField("B",B,zhat,name1));
    builder.addInteraction(interactions.getMagneticField("B",B,zhat,name3));

    SpinWave test = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(1.0,0.0,0.0);
    Line.setFinalPoint(2.0,0.0,0.0);
    Line.setNumberPoints(11);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    SpinWaveDispersion dispersion;
    dispersion.setFilename("VillainFMPhase.txt");
    dispersion.setGenie(test);
    dispersion.setPoints(kPoints);
    
    dispersion.save();
    
    std::cout << "Analytical Results: For comparison with VillainFMPhase.txt\n";
    std::cout << "H, K, L, E1, E+, E-, I+, I-\n";

    double ky = 0.0;
    for(const auto& point : kPoints)
    {
      cout << point.get<0>() << " " << point.get<1>() << " " << point.get<2>() << " ";
      double kx = 2.0*M_PI*point.get<0>();
      double ky = 2.0*M_PI*point.get<0>();
      double R1K = sqrt(pow(eta+1,2)*pow(1-cos(kx),2)+4.0*pow(gamma*cos(ky),2));
      double term = B + 2.0*gamma + (eta-1.0)*(cos(kx)-1.0);
      cout << term+R1K << " " << term-R1K << " ";
      cout << (R1K-2.0*gamma*cos(ky))/R1K << " " << (R1K+2.0*gamma*cos(ky))/R1K << endl;
    }

    return 0;
}
