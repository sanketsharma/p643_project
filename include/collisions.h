#ifndef _COLLISIONS_H
#define _COLLISIONS_H

#include <vector>
#include <Cell.h>
#include <CollisionPartnersGenerator.h>
#include <map>
#include <PostCollisionVelocitiesGenerator.h>
#include <interpolation.h>

namespace p643
{
    double  collisionDepletion(const std::vector<std::vector<std::vector<double>>>& distributionFunctionGrid, 
                                                                const std::vector<Index>& mj, 
                                                                const Index& etaIHat, 
                                                                const Index& zetaIHat, 
                                                                double beta,
                                                                double deltaT);
    void processCollisions(Cell& cell, const std::vector<Index> mj, 
                            PostCollisionVelocitiesGenerator& postCollisionVelocitiesGenerator, 
                            const std::array<double, 3> eta, 
                            const Index& etaIHat, 
                            const Index& zetaIHat, 
                            double beta, 
                            double deltaT);
}

#endif
