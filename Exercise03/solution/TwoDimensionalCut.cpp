#include "CommonFunctions.h"

using namespace std;
using namespace SpinWaveGenie;

int main()
{
    SpinWave SW = createModel();

    PointsAlongLine Line;
    Line.setFirstPoint(3.0,0.0,-4.0);
    Line.setFinalPoint(3.0,0.0, 4.0);
    Line.setNumberPoints(101);
    ThreeVectors<double> kPoints = Line.getPoints();

    Energies energies(0.0, 80.0, 101);

    TwoDimGaussian resinfo;
    resinfo.a = 579.7;
    resinfo.b = -20.0;
    resinfo.c = 1.3;
    resinfo.tol = 1.0e-1;
    resinfo.direction = Vector3(0.0,0.0,1.0);

    unique_ptr<SpinWavePlot> res(new TwoDimensionResolutionFunction(resinfo, SW, energies));

    HKLDirections direction;
    direction.addDirection(1.0,0.0,0.0,0.2);
    direction.addDirection(0.0,1.0,0.0,0.2);
    //direction.addDirection(0.0,0.0,1.0,0.05);
    
    unique_ptr<SpinWavePlot> asdf(new IntegrateAxes(std::move(res),direction,1.0e-1));
   
    TwoDimensionalCut twodimcut;
    twodimcut.setFilename("YFeO3_3_0_L");
    twodimcut.setPlotObject(move(asdf));
    twodimcut.setPoints(kPoints);
    twodimcut.save();
   
    Line.setFirstPoint(3.0,1.0,-4.0);
    Line.setFinalPoint(3.0,1.0, 4.0);
    Line.setNumberPoints(101);
    kPoints = Line.getPoints();
    twodimcut.setFilename("YFeO3_3_1_L");
    twodimcut.setPoints(kPoints);
    twodimcut.save();
    return 0;
}
