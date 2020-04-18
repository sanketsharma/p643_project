#include <Collider.h>
#include <easylogging++.h>

namespace p643
{
Collider::Collider(Cell& cell, double beta, double deltaT):
    cell(cell),
    postCollisionVelocitiesGenerator(),
    beta(beta),
    deltaT(deltaT),
    velocityGridSide(cell.distributionFunctionGrid.size())
{}

double  Collider::collisionDepletion( const std::vector<Index>& mj, 
                                        const Index& etaIHat, 
                                        const Index& zetaIHat)
{
    const std::vector<std::vector<std::vector<double>>>& distributionFunctionGrid = cell.distributionFunctionGrid;
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
    const double phiHatZetaIHat = distributionFunctionGrid[i][j][k];
    
    return 1.0/(2 * mj.size()) * deltaT * beta * beta * beta * sum * sgn(phiHatZetaIHat);          
}

void Collider::replenish(const std::map<Index, double>& replenishmentFractions, double replenishment)
{
    for(const auto& pair : replenishmentFractions)
    {
        const auto& index = pair.first;    
        double fraction = pair.second;    
        int i = std::get<0>(index) + velocityGridSide/2;
        int j = std::get<1>(index) + velocityGridSide/2;
        int k = std::get<2>(index) + velocityGridSide/2;
        cell.distributionFunctionGrid[i][j][k] += fraction*replenishment;
    }
}

void Collider::processCollisions(const std::vector<Index> mj, 
                        const std::array<double, 3> eta, 
                        const Index& etaIHat, 
                        const Index& zetaIHat)
{
    const double depletion = collisionDepletion(mj, etaIHat, zetaIHat);
    const double zetaX = cell.getVelocity(beta, std::get<0>(zetaIHat));
    const double zetaY = cell.getVelocity(beta, std::get<1>(zetaIHat));
    const double zetaZ = cell.getVelocity(beta, std::get<2>(zetaIHat));
    const std::array<double, 3> zeta{{zetaX, zetaY, zetaZ}};

    const auto postCollisionVelocities = postCollisionVelocitiesGenerator.getPostCollisionVelocities(eta, zeta);
    const double etaPrimeX = postCollisionVelocities[0];
    const double etaPrimeY = postCollisionVelocities[1];
    const double etaPrimeZ = postCollisionVelocities[2];
    auto replishmentFractionsEta = interpolateToGrid(velocityGridSide, beta, etaPrimeX, etaPrimeY, etaPrimeZ);
    const double replenishmentEta = depletion/2;
    replenish(replishmentFractionsEta, replenishmentEta); 

    const double zetaPrimeX = postCollisionVelocities[3];
    const double zetaPrimeY = postCollisionVelocities[4];
    const double zetaPrimeZ = postCollisionVelocities[5];
    auto replishmentFractionsZeta = interpolateToGrid(velocityGridSide, beta, zetaPrimeX, zetaPrimeY, zetaPrimeZ);
    const double replenishmentZeta = depletion/2;
    replenish(replishmentFractionsZeta, replenishmentZeta);
}
}
