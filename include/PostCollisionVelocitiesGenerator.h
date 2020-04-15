#ifndef _POSTCOLLISIOnVELOCITIESGENERATOR_H_
#define _POSTCOLLISIOnVELOCITIESGENERATOR_H_

#include "constants.h"
#include <array>
#include <random>

namespace p643
{
    class PostCollisionVelocitiesGenerator
    {
        private:
        std::random_device rd;
        std::default_random_engine e;
        std::uniform_real_distribution<> dis;

        std::array<double, 3> getRelativeVelocity(const std::array<double, 3>& v1, const std::array<double, 3>& v2);
        double getLength(const std::array<double, 3>& v);

        public:
        PostCollisionVelocitiesGenerator();
        std::array<double, 6> getPostCollisionVelocities(const std::array<double, 3>& v1, const std::array<double, 3>& v2);
    };
}

#endif
