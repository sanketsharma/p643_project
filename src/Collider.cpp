#include <Collider.h>
#include <easylogging++.h>
#include <cmath>

namespace p643
{
Collider::Collider(Cell& cell, double beta, double deltaT):
    cell(cell),
    postCollisionVelocitiesGenerator(),
    beta(beta),
    deltaT(deltaT),
    velocityGridSide(cell.distributionFunctionGrid.size())
{}

Collider::~Collider()
{
}

double  Collider::collisionDepletion( const std::vector<Index>& mj, 
                                        const Index& etaIHat)
{
    const std::vector<std::vector<std::vector<double>>>& distributionFunctionGrid = cell.distributionFunctionGrid;
    double sum = 0.0;
    double phiHatEtaIHat = 0.0;
    for(const auto& index : mj)
    {
        const int x = std::get<0>(index);
        const int y = std::get<1>(index);
        const int z = std::get<2>(index);

        if(index == etaIHat)
        {
            phiHatEtaIHat = std::abs(distributionFunctionGrid[x][y][z]);
        }
        sum += std::abs(distributionFunctionGrid[x][y][z]);
    }

    const int p = std::get<0>(etaIHat);
    const int q = std::get<1>(etaIHat);
    const int r = std::get<2>(etaIHat);
    double retVal = 1.0/(2 * sum) * distributionFunctionGrid[p][q][r] * deltaT * beta * beta * beta * (sum - phiHatEtaIHat);
    return retVal;          
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

void Collider::processCollisions( const std::array<double, 3> eta, 
                        const Index& etaIHat, 
                        const Index& zetaIHat, 
                        const double depletionAbsValue)
{
    const int zetaIHatX = std::get<0>(zetaIHat);
    const int zetaIHatY = std::get<1>(zetaIHat);
    const int zetaIHatZ = std::get<2>(zetaIHat);
    const double phiHatZetaIHat = cell.distributionFunctionGrid[zetaIHatX][zetaIHatY][zetaIHatZ];
    const double depletion = depletionAbsValue * sgn(phiHatZetaIHat);

    const double zetaX = cell.getVelocity(beta, zetaIHatX);
    const double zetaY = cell.getVelocity(beta, zetaIHatY);
    const double zetaZ = cell.getVelocity(beta, zetaIHatZ);
    const std::array<double, 3> zeta{{zetaX, zetaY, zetaZ}};

    const auto postCollisionVelocities = postCollisionVelocitiesGenerator.getPostCollisionVelocities(eta, zeta);
    const double etaPrimeX = postCollisionVelocities[0];
    const double etaPrimeY = postCollisionVelocities[1];
    const double etaPrimeZ = postCollisionVelocities[2];
    auto replenishmentFractionsEta = interpolateToGrid(velocityGridSide, beta, etaPrimeX, etaPrimeY, etaPrimeZ);
    if(!replenishmentFractionsEta.empty()) //In case no replenishment is done we won't do depletion as well
    {
        int etaIHatX = std::get<0>(etaIHat);
        int etaIHatY = std::get<1>(etaIHat);
        int etaIHatZ = std::get<2>(etaIHat);
        cell.distributionFunctionGrid[etaIHatX][etaIHatY][etaIHatZ] -= depletion/2;
        const double replenishmentEta = depletion/2;
        replenish(replenishmentFractionsEta, replenishmentEta); 
    }
    const double zetaPrimeX = postCollisionVelocities[3];
    const double zetaPrimeY = postCollisionVelocities[4];
    const double zetaPrimeZ = postCollisionVelocities[5];
    auto replenishmentFractionsZeta = interpolateToGrid(velocityGridSide, beta, zetaPrimeX, zetaPrimeY, zetaPrimeZ);
    if(!replenishmentFractionsZeta.empty())//In case no replenishment is done we won't do depletion as well
    {
        cell.distributionFunctionGrid[zetaIHatX][zetaIHatY][zetaIHatZ] -= depletion/2;
        const double replenishmentZeta = depletion/2;
        replenish(replenishmentFractionsZeta, replenishmentZeta);
    }
}
}
