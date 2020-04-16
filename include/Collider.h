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
class Collider
{
    public:
    Collider(Cell& cell, double beta, double deltaT);
    void processCollisions(const std::vector<Index> mj, 
                            const std::array<double, 3> eta, 
                            const Index& etaIHat, 
                            const Index& zetaIHat);

    private:
    double  collisionDepletion( const std::vector<Index>& mj, 
                                const Index& etaIHat, 
                                const Index& zetaIHat);
    Cell& cell;
    PostCollisionVelocitiesGenerator postCollisionVelocitiesGenerator;
    void replenish(const std::map<Index, double>& replenishmentFractions, double replenishment);
    const double beta;
    const double deltaT;
    const int velocityGridSide;

};
}

#endif
