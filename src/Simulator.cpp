#include <Collider.h>
#include <cmath>
#include <tuple>
#include <CollisionPartnersGenerator.h>
#include <Simulator.h>
#include <interpolation.h>
#include <easylogging++.h>

namespace p643
{

Simulator::Simulator(const std::map<std::string, std::string>& configValues):
    k(std::stod(configValues.at("k"))),
    Tr(std::stod(configValues.at("Tr"))),
    m(std::stod(configValues.at("m"))),
    etaR(std::sqrt(2*k*Tr/m)),
    u(std::sqrt(8*k*Tr/(p643::pi*m))),
    n(std::stoi(configValues.at("n"))),
    size(std::stod(configValues.at("size"))),
    beta(std::stod(configValues.at("beta"))),
    tHat(1.0),
    uHat(u/etaR),
    gridX(std::stoi(configValues.at("gridX"))),
    gridY(std::stoi(configValues.at("gridY"))),
    gridZ(std::stoi(configValues.at("gridZ"))),
    depletingFraction(std::stod(configValues.at("depletingFraction"))),
    velocityGridSide(std::stoi(configValues.at("velocityGridSide"))),
    maximumIndex(velocityGridSide),
    grid(gridX, gridY, gridZ, velocityGridSide, size, beta, etaR, n, tHat, uHat),
    postCollisionVelocitiesGenerator(),
    deltaT(std::stod(configValues.at("deltaT"))),
    maxSteps(100)
{
}

void Simulator::simulate()
{
    for(int step = 0 ; step < maxSteps; ++step)
    {
        for(int p = 0; p < gridX ; ++p)
        {
            for(int q = 0; q < gridY; ++q)
            {
                for(int r = 0; r < gridZ; ++r)
                {
                    auto& cell = grid.myGrid[p][q][r];
                    CollisionPartnersGenerator collisionPartnersGenerator(cell);
                    Collider collider(cell, beta, deltaT);

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
                                    collider.processCollisions(mj, eta, etaIHat, zetaIHat);
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
