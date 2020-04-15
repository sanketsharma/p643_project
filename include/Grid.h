#ifndef _GRID_H_
#define _GRID_H_

#include <Cell.h>
#include <constants.h>
#include <vector>

namespace p643
{
//Division of space
struct Grid
{
    Grid(const unsigned l, const unsigned m, const unsigned n, const unsigned velocityGridSide, const double size, const double beta, const double etaR, const double nHat, const double tHat, const double uHat);
    std::vector<std::vector<std::vector<Cell>>> myGrid;
};
}

#endif
