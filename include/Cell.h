#ifndef CELL_H 
#define CELL_H 

#include <constants.h>
#include <vector>

namespace p643
{
    struct Cell
    {
        const unsigned velocityGridSide;
        const double size;
        const double beta;
        const double etaR;
        const int nHat;
        const double tHat;
        const double uHat;
        std::vector<std::vector<std::vector<double>>> distributionFunctionGrid;

        Cell(const unsigned velocityGridSide, const double size, const double beta, const double etaR, const double nHat, const double tHat, const double uHat);

        static inline double getVelocity(const double beta, const unsigned index)
        {
            return beta * (index - 6);
        }
    };
}

#endif
