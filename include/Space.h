#ifndef _SPACE_H_
#define _SPACE_H_

#include <array>
#include <Cell.h>
namespace p643
{
class Space
{
    const int N;
    std::array<std::array<Cell, N>, N> Cells;

    public:
    Space(int N): N(N)
    {}
};

}

#endif
