#include <boost/test/unit_test.hpp>
#include <interpolation.h>
#include <iostream>

double testHelper(int velocityGridSide, double beta, double vx, double vy, double vz)
{
    std::map<p643::Index, double> interpolatedResults = p643::interpolateToGrid(velocityGridSide/2, beta, vx, vy, vz);
    double sum = 0;
    for(auto e: interpolatedResults)
    {
        auto f = p643::toVelocity(e.first, velocityGridSide);
        int x = std::get<0>(f);
        int y = std::get<1>(f);
        int z = std::get<2>(f);
        double w   = e.second;
        sum += w;
        std::string delim = ", " ;
        std::cout << x ;
        std::cout << delim; 
        std::cout << y;
        std::cout << delim;
        std::cout<< z ;
        std::cout << delim;
        std::cout << w;
        std::cout << std::endl;
    }
    return sum;
}
BOOST_AUTO_TEST_CASE(interpolation_test)
{

    const unsigned velocityGridSide = 11; 
    const double beta = 0.6; 

    double vx = 4.1/beta;
    double vy = 4.1/beta;
    double vz = 4.1/beta;
    int sum = testHelper(velocityGridSide, beta, vx, vy, vz);
    BOOST_CHECK_CLOSE(sum, 1.0, 0.001);

    std::cout << std::endl;

    vx = 5.1/beta;
    vy = 4.1/beta;
    vz = 4.1/beta;
    sum = testHelper(velocityGridSide, beta, vx, vy, vz);
    BOOST_CHECK_CLOSE(sum, 1.0, 0.001);

    std::cout << std::endl;

    vx = 5.1/beta;
    vy = 5.1/beta;
    vz = 5.1/beta;
    testHelper(velocityGridSide, beta, vx, vy, vz);
    sum = testHelper(velocityGridSide, beta, vx, vy, vz);
    BOOST_CHECK_CLOSE(sum, 0.0, 0.001);

}
