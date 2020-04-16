#ifndef __DISTRIBUTION__H
#define __DISTRIBUTION__H

#include <distribution.h>
#include <cmath>
namespace p643
{
    double getEquilibriumDistribution(double etaHat, double nHat, double uHat, double tHat)
    {
        return nHat * std::pow(pi*tHat, -3.0/2.0) * std::exp(-1.0/tHat*std::pow(etaHat-uHat, 2.0));
    }
}

#endif
