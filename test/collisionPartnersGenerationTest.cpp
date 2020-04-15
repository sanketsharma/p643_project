#define BOOST_TEST_MODULE DVMTests

#include <boost/test/unit_test.hpp>
#include <Grid.h>
#include <CollisionPartnersGenerator.h>

BOOST_AUTO_TEST_CASE(collision_partners_generator_test)
{
    const double k = 1.38064852e-23;
    const double Tr = 300;   //Kelvin
    const double m = 6.6335209e-26; //kg

    const double etaR = std::sqrt(2*k*Tr/m);
    const double u = std::sqrt(8*k*Tr/(p643::pi*m));
    const int n = 10000000;
    const double size = 1e-6; //meters
    const double beta = 0.6;

    double tHat = 1;
    double uHat = u/etaR;

 
    p643::Grid grid(1, 1, 1, 11, size, beta, etaR, n, tHat, uHat);
    p643::Cell& cell = grid.myGrid[0][0][0];
    p643::CollisionPartnersGenerator collisionPartnersGenerator(cell);
    double depletionFactor = 0.0;
    auto shouldContainOne = collisionPartnersGenerator.getCollisionPartners(depletionFactor);
    BOOST_CHECK(shouldContainOne.size() == 1) ;

    depletionFactor = 0.5;
    auto shouldContainAboutHalf = collisionPartnersGenerator.getCollisionPartners(depletionFactor);
}
