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
