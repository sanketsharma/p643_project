#include <CollisionPartnersGenerator.h>
#include <algorithm>
#include <easylogging++.h>

namespace p643
{

CollisionPartnersGenerator::CollisionPartnersGenerator(const Cell& cell):
distributionFunctionGrid(cell.distributionFunctionGrid),
sum(0.0),
i(cell.distributionFunctionGrid.size()),
j(cell.distributionFunctionGrid.begin()->size()),
k(cell.distributionFunctionGrid.begin()->begin()->size()),
nV(i*j*k),
cumulativeDistributionFunction(nV),
generator(randomDevice())
{}

void CollisionPartnersGenerator::populateCumulativeDistributionFunction()
{
    auto it = cumulativeDistributionFunction.begin();
    sum = 0.0;
    for(auto it1 = distributionFunctionGrid.begin(); it1 != distributionFunctionGrid.end(); ++it1)
    {
        for(auto it2 = it1->begin(); it2 != it1->end(); ++it2)
        {
            for(auto it3 = it2->begin(); it3 != it2->end(); ++it3)
            {
                sum += std::abs(*it3);
                *it = sum;
                ++it;
            }
        }
    }
    decltype(distribution.param()) newRange(0.0, sum);
    distribution.param(newRange);
}

std::vector<Index> CollisionPartnersGenerator::getCollisionPartners(double depletingFraction)
{
    double cumulativeMassFraction = 0.0;
    std::vector<Index> collisionPartners;
    populateCumulativeDistributionFunction();

    while(true)
    {
        const double randomNumber = distribution(generator);
        auto it = std::lower_bound(cumulativeDistributionFunction.begin(), cumulativeDistributionFunction.end(), randomNumber);

        if(it != cumulativeDistributionFunction.begin())
        {
            --it;
        }
              
        const unsigned l = it - cumulativeDistributionFunction.begin();
        const unsigned x = l/(i*j);
        const unsigned y = (l%(i*j))/j;
        const unsigned z = (l%(i*j))%j;

        collisionPartners.emplace_back(x,y,z);

        cumulativeMassFraction += std::abs(distributionFunctionGrid[x][y][z]/sum);

        if(cumulativeMassFraction >= depletingFraction || collisionPartners.size() >= nV) 
        {
            break;
        }
    }
    return collisionPartners;
}

}
