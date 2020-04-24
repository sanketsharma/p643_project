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
    void processCollisions( const std::array<double, 3> eta, 
                            const Index& etaIHat, 
                            const Index& zetaIHat, 
                            const double depletionAbsValue);
    ~Collider();

    double  collisionDepletion( const std::vector<Index>& mj, 
                                const Index& etaIHat);
    private:
    void replenish(const std::map<Index, double>& replenishmentFractions, double replenishment);

    Cell& cell;
    PostCollisionVelocitiesGenerator postCollisionVelocitiesGenerator;
    const double beta;
    const double deltaT;
    const int velocityGridSide;
};
}

#endif
