#include <Collider.h>
#include <cmath>
#include <tuple>
#include <chrono>
#include <fstream>
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
    grid(gridX, gridY, gridZ, velocityGridSide, size, beta, etaR, n, tHat, uHat),
    postCollisionVelocitiesGenerator(),
    deltaT(std::stod(configValues.at("deltaT"))),
    maxSteps(std::stoi(configValues.at("maxSteps")))
{}

void Simulator::dumpData()
{
    std::map<int, double> flattenedDistribution;
    for(int p = 0; p < gridX; ++p)
    {
        for(int q = 0; q < gridY; ++q)
        {
            for(int r = 0; r < gridZ; ++r)
            {
                const auto& cell = grid.myGrid[p][q][r];
                const auto& distributionFunctionGrid = cell.distributionFunctionGrid;
                int shift = velocityGridSide/2;
                for(int i = 0; i < velocityGridSide; ++i)
                {
                    for(int j = 0; j < velocityGridSide; ++j)
                    {
                        for(int k = 0; k < velocityGridSide; ++k)
                        {
                            const int indexSquaredSum = (i-shift)*(i-shift) + (j-shift)*(j-shift) + (k-shift)*(k-shift);
                            auto it = flattenedDistribution.find(indexSquaredSum);
                            if(it != flattenedDistribution.end())
                            {
                                flattenedDistribution[indexSquaredSum] += distributionFunctionGrid[i][j][k];
                            }
                            else
                            {
                                flattenedDistribution.emplace(indexSquaredSum, distributionFunctionGrid[i][j][k]);
                            }
                        }
                    }
                }
            }
        }
    }

    long ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    std::ofstream ofs(std::string("../logs/cell.dat.")+ std::to_string(ms));
    for(auto pair : flattenedDistribution)
    {
        ofs << beta * std::sqrt(pair.first)*etaR << " " << pair.second << std::endl;
    }
    ofs.close();
}

void Simulator::simulate()
{
    dumpData();
    for(int step = 0; step < maxSteps; ++step)
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

                    for(int i = 0; i < velocityGridSide; ++i)
                    {
                        const double etaX = cell.getVelocity(beta, i);
                        for(int j = 0; j < velocityGridSide; ++j)
                        {
                            const double etaY = cell.getVelocity(beta, j);
                            for(int k = 0; k < velocityGridSide; ++k)
                            {
                                const double etaZ = cell.getVelocity(beta, k);
                                const auto etaIHat = std::make_tuple(i, j, k);
                                const auto mj = collisionPartnersGenerator.getCollisionPartners(depletingFraction);
                                if(mj.empty())
                                {
                                    continue;
                                }
                                const std::array<double, 3> eta{{etaX, etaY, etaZ}};
                                const double depletionAbsValue = collider.collisionDepletion(mj, etaIHat);

                                for(auto zetaIHat :  mj)
                                {
                                    collider.processCollisions(eta, etaIHat, zetaIHat, depletionAbsValue);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    dumpData();
}
}
