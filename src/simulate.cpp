#include <Grid.h>
#include <Cell.h>
#include <collisions.h>
#include <cmath>
#include <tuple>
#include <CollisionPartnersGenerator.h>
#include <PostCollisionVelocitiesGenerator.h>
#include <simulate.h>
#include <interpolation.h>
#include <easylogging++.h>

namespace p643
{

void replenish(Cell& cell, const std::map<Index, double>& replenishmentFractions, double replenishment)
{
    for(const auto& pair : replenishmentFractions)
    {
        const auto& index = pair.first;    
        double fraction = pair.second;    
        int i = std::get<0>(index);
        int j = std::get<1>(index);
        int k = std::get<2>(index);
        cell.distributionFunctionGrid[i][j][k] += fraction*replenishment;
    }
}

void simulate(const std::map<std::string, std::string>& configValues)
{
    const double k = std::stod(configValues.at("k"));
    const double Tr = std::stod(configValues.at("Tr"));
    const double m = std::stod(configValues.at("m"));
    const double etaR = std::sqrt(2*k*Tr/m);
    const double u = std::sqrt(8*k*Tr/(p643::pi*m));
    const int n = std::stoi(configValues.at("n"));
    const double size = std::stod(configValues.at("size"));
    const double beta = std::stod(configValues.at("beta"));
    const double tHat = 1.0;
    const double uHat = u/etaR;
    const int gridX = std::stoi(configValues.at("gridX"));
    const int gridY = std::stoi(configValues.at("gridY"));
    const int gridZ = std::stoi(configValues.at("gridZ"));
    const double depletingFraction = std::stod(configValues.at("depletingFraction"));
    const unsigned velocityGridSide = std::stoi(configValues.at("velocityGridSide"));
    const unsigned maximumIndex = 2*velocityGridSide + 1;

    p643::Grid grid(gridX, gridY, gridZ, velocityGridSide, size, beta, etaR, n, tHat, uHat);
    p643::PostCollisionVelocitiesGenerator postCollisionVelocitiesGenerator;

    double deltaT = std::stod(configValues.at("deltaT"));
    for(unsigned step = 0 ; step < 100; ++step)
    {
        for(int p = 0; p < gridX ; ++p)
        {
            for(int q = 0; q < gridY; ++q)
            {
                for(int r = 0; r < gridZ; ++r)
                {
                    auto& cell = grid.myGrid[p][q][r];
                    p643::CollisionPartnersGenerator collisionPartnersGenerator(cell);

                    for(unsigned i = 0; i < maximumIndex; ++i)
                    {
                        const double etaX = p643::Cell::getVelocity(beta, i);
                        for(unsigned j = 0; j < maximumIndex; ++j)
                        {
                            const double etaY = p643::Cell::getVelocity(beta, j);
                            for(unsigned k = 0; k < maximumIndex; ++k)
                            {
                                const double etaZ = p643::Cell::getVelocity(beta, k);
                                const auto etaIHat = std::make_tuple(i, j, k);
                                const auto mj = collisionPartnersGenerator.getCollisionPartners(depletingFraction);
                                const std::array<double, 3> eta{{etaX, etaY, etaZ}};

                                for(auto zetaIHat :  mj)
                                {
                                    const double depletion = p643::collisionDepletion(cell.distributionFunctionGrid, mj, etaIHat, zetaIHat, beta, deltaT);
                                    const double zetaX = p643::Cell::getVelocity(beta, std::get<0>(zetaIHat));
                                    const double zetaY = p643::Cell::getVelocity(beta, std::get<1>(zetaIHat));
                                    const double zetaZ = p643::Cell::getVelocity(beta, std::get<2>(zetaIHat));
                                    const std::array<double, 3> zeta{{zetaX, zetaY, zetaZ}};

                                    const auto postCollisionVelocities = postCollisionVelocitiesGenerator.getPostCollisionVelocities(eta, zeta);
                                    const double etaPrimeX = postCollisionVelocities[0];
                                    const double etaPrimeY = postCollisionVelocities[1];
                                    const double etaPrimeZ = postCollisionVelocities[2];
                                    auto replishmentFractionsEta = interpolateToGrid(velocityGridSide, beta, etaPrimeX, etaPrimeY, etaPrimeZ);
                                    const double replenishmentEta = depletion/2;
                                    replenish(cell, replishmentFractionsEta, replenishmentEta);

                                    const double zetaPrimeX = postCollisionVelocities[3];
                                    const double zetaPrimeY = postCollisionVelocities[4];
                                    const double zetaPrimeZ = postCollisionVelocities[5];
                                    auto replishmentFractionsZeta = interpolateToGrid(velocityGridSide, beta, zetaPrimeX, zetaPrimeY, zetaPrimeZ);
                                    const double replenishmentZeta = depletion/2;
                                    replenish(cell, replishmentFractionsZeta, replenishmentZeta);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
}
