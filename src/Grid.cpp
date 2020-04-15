#include <Grid.h>

namespace p643
{
    Grid::Grid(const unsigned l, const unsigned m, const unsigned n, const unsigned velocityGridSide, const double size, const double beta, const double etaR, const double nHat, const double tHat, const double uHat):
    myGrid(l, std::vector<std::vector<Cell>>(m, std::vector<Cell>(n, Cell(velocityGridSide, size, beta, etaR, nHat, tHat, uHat))))
    {}
}

