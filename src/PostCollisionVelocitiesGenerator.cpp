#include <PostCollisionVelocitiesGenerator.h>
#include <cmath>
#include <numeric>

namespace p643
{
    
    PostCollisionVelocitiesGenerator::PostCollisionVelocitiesGenerator()
    :e(rd()),
    dis(0.0, 1.0)
    {}
    std::array<double, 3> PostCollisionVelocitiesGenerator::getRelativeVelocity(const std::array<double, 3>& v1, 
                                                                                const std::array<double, 3>& v2)
    {
        return std::array<double, 3> {{v1[0]-v2[0], v1[1]-v1[1], v1[2]-v2[2]}};
    }

    double PostCollisionVelocitiesGenerator::getLength(const std::array<double, 3>& v) 
    {
        double x = std::accumulate(v.begin(), v.end(), 0.0, [](double a, double b){
                                                                          return a + std::pow(b, 2.0);
                                                                       });
        return std::sqrt(x);
    }

    //See https://mathworld.wolfram.com/SpherePointPicking.html
    std::array<double, 6> PostCollisionVelocitiesGenerator::getPostCollisionVelocities(const std::array<double, 3>& eta, 
                                                                                        const std::array<double, 3>& zeta)
    {
        double u = dis(e);
        double theta = 2*pi*u;

        double v = dis(e);
        double phi = std::acos(2*v - 1);

        const std::array<double, 3> relativeVelocity =  getRelativeVelocity(eta, zeta);
        const std::array<double, 3> cMVelocity{{(eta[0]+zeta[0])/2.0, (eta[1]+zeta[1])/2.0, (eta[2]+zeta[2])/2.0}};
        double radius = getLength(relativeVelocity)/2.0;
        const std::array<double, 3> temp ={{ radius*std::cos(theta)*std::sin(phi), radius*std::sin(theta)*std::sin(phi), 
                                                                                    radius*std::cos(phi)}};

        const std::array<double, 3> etaPrime = {{cMVelocity[0]+temp[0], cMVelocity[1]+temp[1], cMVelocity[2]+temp[2]}};
        const std::array<double, 3> zetaPrime = {{cMVelocity[0]-temp[0], cMVelocity[1]-temp[1], cMVelocity[2]-temp[2]}};

        return std::array<double, 6> {{etaPrime[0], etaPrime[1], etaPrime[2], zetaPrime[0], zetaPrime[1], zetaPrime[2]}};
    }
}
