#include <boost/test/unit_test.hpp>
#include <interpolation.h>

void testHelper(int maxIndexMagnitude, double beta, double vx, double vy, double vz)
{
    std::map<p643::Index, double> interpolatedResults = p643::interpolateToGrid(maxIndexMagnitude, beta, vx, vy, vz);
    for(auto e: interpolatedResults)
    {
        int x = std::get<0>(e.first);
        int y = std::get<1>(e.first);
        int z = std::get<2>(e.first);
        double w   = e.second;
        std::string delim = " , " ;
        std::cout << x ;
        std::cout << delim; 
        std::cout << y;
        std::cout << delim;
        std::cout<< z ;
        std::cout << delim;
        std::cout << w;
        std::cout << std::endl;
    }
}
BOOST_AUTO_TEST_CASE(interpolation_test)
{

    const unsigned maxIndexMagnitude = 5; 
    const double beta = 0.6; 

    double vx = 4.1*beta;
    double vy = 4.1*beta;
    double vz = 4.1*beta;
    testHelper(maxIndexMagnitude, beta, vx, vy, vz);

     std::cout << std::endl;

     vx = 5.1*beta;
     vy = 4.1*beta;
     vz = 4.1*beta;
     testHelper(maxIndexMagnitude, beta, vx, vy, vz);

     std::cout << std::endl;

     vx = 5.1*beta;
     vy = 5.1*beta;
     vz = 5.1*beta;
     testHelper(maxIndexMagnitude, beta, vx, vy, vz);
    //BOOST_CHECK(shouldContainOne.size() == 1) ;

}
