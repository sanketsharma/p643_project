#ifndef _COLLISIONPARTNERS_
#define _COLLISIONPARTNERS_

#include <constants.h>
#include <Cell.h>
#include <vector>
#include <tuple>
#include <random>

namespace p643
{
class CollisionPartnersGenerator
{
    public:
    CollisionPartnersGenerator(const Cell& cell);
    std::vector<Index> getCollisionPartners(double depletingFraction);

    private:
    void populateCumulativeDistributionFunction();

    const std::vector<std::vector<std::vector<double>>>& distributionFunctionGrid;
    double sum;
    const unsigned i;
    const unsigned j;
    const unsigned k;
    const unsigned nV;
    std::vector<double> cumulativeDistributionFunction;
    std::random_device randomDevice;
    std::mt19937 generator;
    std::uniform_real_distribution<> distribution;
};
}

#endif
