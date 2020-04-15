#include <collisions.h>

namespace p643
{
    double  collisionDepletion(const std::vector<std::vector<std::vector<double>>>& distributionFunctionGrid, 
                                                                const std::vector<Index>& mj, 
                                                                const Index& etaIHat, 
                                                                const Index& zetaIHat, 
                                                                double beta,
                                                                double deltaT)
    {
        double sum = 0.0;
        for(const auto& index : mj)
        {
            if(index == etaIHat)
            {
                continue;
            }

            const unsigned x = std::get<0>(index);
            const unsigned y = std::get<1>(index);
            const unsigned z = std::get<2>(index);
            sum += std::abs(distributionFunctionGrid[x][y][z]);
        }

        const unsigned i = std::get<0>(zetaIHat);
        const unsigned j = std::get<1>(zetaIHat);
        const unsigned k = std::get<2>(zetaIHat);
        const double phiHatZetaIHat = distributionFunctionGrid[i][j][j];
        
        return 1.0/(2 * mj.size()) * deltaT * beta * beta * beta * sum * sgn(phiHatZetaIHat);          
    }
}
