#include <Cell.h>
#include <cmath>
#include <distribution.h>
#include <stdexcept>
#include <easylogging++.h>

namespace p643
{

    Cell::Cell(const unsigned velocityGridSide, const double size, const double beta, const double etaR, const double nHat, const double tHat, const double uHat):
    velocityGridSide(velocityGridSide),
    size(size),
    beta(beta),
    etaR(etaR),
    nHat(nHat),
    tHat(tHat),
    uHat(uHat),
    distributionFunctionGrid(velocityGridSide, std::vector<std::vector<double>>(velocityGridSide, std::vector<double>(velocityGridSide, 0.0)))
    {
        if(velocityGridSide%2==0)
        {
            throw std::invalid_argument("Velocity Grid Side can only have odd values");
        }

        for(unsigned i = 0; i < velocityGridSide; ++i)
        {
            double vxHat = getVelocity(beta, i);
            for(unsigned j = 0; j < velocityGridSide; ++j)    
            {
                double vyHat = getVelocity(beta, j);
                for(unsigned k = 0 ; k < velocityGridSide; ++k)
                {
                    double vzHat = getVelocity(beta, k);
                    double speedHat = std::sqrt(vxHat*vxHat + vyHat*vyHat + vzHat*vzHat);
                    distributionFunctionGrid[i][j][k] = getEquilibriumDistribution(speedHat, nHat, uHat, tHat);
                }
            }
        }
    }
}
