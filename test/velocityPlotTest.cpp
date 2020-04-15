#include <PostCollisionVelocitiesGenerator.h>
#include <iostream>
#include <fstream>

int main()
{
    int n = 0;
    std::cout << "Enter number of velocities to generate: ";
    std::cin >> n;
    std::array<double, 3> eta = {{1.0, 1.0, 1.0}};
    std::array<double, 3> zeta = {{2.0, 2.0, 2.0}};

    p643::PostCollisionVelocitiesGenerator g;

    std::ofstream ofs("velocities.dat");
    for(int i =0; i < n; ++i)
    {
        std::array<double,6> postCollisionVelocities = g.getPostCollisionVelocities(eta, zeta);
        std::array<double,3> etaPrime = {{postCollisionVelocities[0], postCollisionVelocities[1], postCollisionVelocities[2]}};
        std::array<double,3> zetaPrime = {{postCollisionVelocities[3], postCollisionVelocities[4], postCollisionVelocities[5]}};
        
        ofs  << etaPrime[0] << ' ' << etaPrime[1] << ' ' << etaPrime[2] << ' ';
        ofs  << zetaPrime[0] << ' ' << zetaPrime[1] << ' ' << zetaPrime[2] << '\n';
    }

    ofs.close();
    return 0;
}
