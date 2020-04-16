#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <tuple>

namespace p643
{
    const double pi = 3.141592653589793238463;
    using Index=std::tuple<int, int, int>;

    inline int sgn(double x)
    {
        return (x > 0) - (x < 0);
    }

    inline Index toVelocity(const Index& i, int gridSide)
    {
        return std::make_tuple(std::get<0>(i) - gridSide/2, std::get<1>(i) - gridSide/2, std::get<2>(i) - gridSide/2);
    }
}

#endif
