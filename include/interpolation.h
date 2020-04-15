#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H

#include <constants.h>
#include <array>
#include <map>
#include <Cell.h>

namespace p643
{
    std::array<double, 5> 
    getFractionalDensityChanges(double a, double b, double c, bool withinVSD, unsigned whichExternalPoint);

    bool isCorner(const unsigned maxIndexMagnitude, Index point);

    std::map<Index, double> 
    interpolateToGrid(const unsigned maxIndexMagnitude, const double beta, const double vx, const double vy, const double vz); //TODO: Test it
}

#endif
