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
    std::array<double,4> BValues{{0.0,4.0,8.0,10.0}};
    for (auto B:BValues)
    {
      double gamma,eta,J;
      gamma = -2.0;
      eta = -2.0;
      J = 1.0;
    
      double theta = M_PI_2;
      theta = acos(-B/(4.0*gamma));
      if (theta != theta) //check for nan
        theta = 0.0;
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
      Line.setFirstPoint(0.5,0.0,0.0);
      Line.setFinalPoint(0.0,0.0,0.0);
      Line.setNumberPoints(101);
      ThreeVectors<double> kPoints = Line.getPoints();
    
      Line.setFirstPoint(0.0,0.0,0.0);
      Line.setFinalPoint(0.0,1.0,0.0);
      Line.setNumberPoints(101);
      ThreeVectors<double> kPoints2 = Line.getPoints();
    
      SpinWaveDispersion dispersion;
      dispersion.setFilename("SC2Chain_"+std::to_string(static_cast<int>(B))+".txt");
      dispersion.setGenie(test);
      dispersion.setPoints(kPoints);
      dispersion.setPoints(kPoints2);
      dispersion.save();
    
    }
    return 0;
}
