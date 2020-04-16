#ifndef _SIMULATE_H
#define _SIMULATE_H

#include <map>
#include <string>
#include <Grid.h>
#include <Cell.h>
#include <PostCollisionVelocitiesGenerator.h>

namespace p643{
class Simulator
{
    public:
    Simulator(const std::map<std::string, std::string>& configValues);
    void simulate();

    private:
    const double k;
    const double Tr;
    const double m;
    const double etaR;
    const double u;
    const int n;
    const double size;
    const double beta;
    const double tHat;
    const double uHat;
    const int gridX;
    const int gridY;
    const int gridZ;
    const double depletingFraction;
    const unsigned velocityGridSide;
    const unsigned maximumIndex;
    Grid grid;
    PostCollisionVelocitiesGenerator postCollisionVelocitiesGenerator;
    double deltaT;
    int maxSteps;

};
} 

#endif
