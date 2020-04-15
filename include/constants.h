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
}

#endif
