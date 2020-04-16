#include <interpolation.h>
#include <cmath>
#include <numeric>
#include <easylogging++.h>


namespace p643
{
std::array<double, 5> 
getFractionalDensityChanges(double a, double b, double c, bool withinVSD, unsigned whichExternalPoint)
{
    const double aSquared = a*a;
    const double bSquared = b*b;
    const double cSquared = c*c;

    const double fo = 1 - aSquared - bSquared - cSquared;
    if(withinVSD)
    {
        const double fext = -(a + b + c - aSquared - bSquared - cSquared)/6;
        const double fix = a + fext;
        const double fiy = b + fext;
        const double fiz = c + fext;
        return std::array<double, 5>{{fo, fext, fix, fiy, fiz}};
    }
    else
    {
        const double fext = -(a + b + c - aSquared - bSquared - cSquared)/2;
        double fix = a;
        double fiy = b;
        double fiz = c;
        if(whichExternalPoint == 0)
        {
            fix += fext;
        }
        else if(whichExternalPoint == 1)
        {
            fix += fext;
        }
        else if(whichExternalPoint == 2)
        {
            fix += fext;
        }
        return std::array<double, 5>{{fo, fext, fix, fiy, fiz}};
    }
}

bool isCorner(const int maxIndexMagnitude, Index point)
{
    int x = std::get<0>(point);
    int y = std::get<1>(point);
    int z = std::get<2>(point);

    if(std::abs(x) == maxIndexMagnitude - 1)
    {
        if(std::abs(y) == maxIndexMagnitude - 1)
        {
            if(std::abs(z) == maxIndexMagnitude - 1)
            {
                return true;
            }
        }

    }
    return false;
}

std::array<unsigned, 3> pointLocation(const int maxIndexMagnitude, const double vxNormalized, const double vyNormalized, const double vzNormalized)
{
    std::array<unsigned, 3> location = {0};

    if(vxNormalized >= maxIndexMagnitude || vxNormalized <= -1*maxIndexMagnitude) //Only one can happen
    {
        location[0] = 1;
    }
    
    if(vyNormalized >= maxIndexMagnitude || vyNormalized <= -1*maxIndexMagnitude)//Only one can happen
    {
        location[1] = 1;
    }
    
    if(vzNormalized >= maxIndexMagnitude || vzNormalized <= -1*maxIndexMagnitude)//Only one can happen
    {
        location[2] = 1;
    }

    return location; //0 => inside 1=> outside face 2=> outside edge 3=> outside corner
}

std::array<Index, 5> getStencilOutsidePoint(const int maxIndexMagnitude, const double vxNormalized, const double vyNormalized, const double vzNormalized, const std::array<unsigned,3>& location, double& a, double& b, double& c, unsigned& whichExternalPoint)
{
        unsigned locationSum = std::accumulate(location.begin(), location.end(), 0);
        Index origin;
        {
            int originX =  static_cast<int>(std::nearbyint(vxNormalized)) + 5;
            int originY =  static_cast<int>(std::nearbyint(vyNormalized)) + 5;
            int originZ =  static_cast<int>(std::nearbyint(vzNormalized)) + 5;
            origin = std::make_tuple(originX, originY, originZ);
        }
        Index ix = origin;
        Index iy = origin;
        Index iz = origin;
        Index e = origin;
        if(locationSum == 1) //Face
        {
            if(location[0] == 1)
            {
                a = vxNormalized - maxIndexMagnitude;
                std::get<0>(origin) = sgn(a) * (maxIndexMagnitude - 1);
                b = vyNormalized - std::get<1>(origin);
                c = vzNormalized - std::get<2>(origin);
                std::get<0>(ix) = std::get<0>(origin) - sgn(a);
                std::get<1>(iy) = std::get<1>(origin) + sgn(b);
                std::get<2>(iz) = std::get<2>(origin) + sgn(c);
                std::get<1>(e) = std::get<1>(origin) - sgn(b);
                whichExternalPoint = 1;
            }
            else if(location[1] == 1)
            {
                b = vyNormalized - maxIndexMagnitude;
                std::get<1>(origin) = sgn(b) * (maxIndexMagnitude - 1);
                c = vzNormalized - std::get<2>(origin);
                a = vxNormalized - std::get<0>(origin);
                std::get<0>(ix) = std::get<0>(origin) + sgn(a);
                std::get<1>(iy) = std::get<1>(origin) - sgn(b);
                std::get<2>(iz) = std::get<2>(origin) + sgn(c);
                std::get<2>(e) = std::get<2>(origin) - sgn(c);
                whichExternalPoint = 2;
            }
            else //location[2] == 1
            {
                c = vzNormalized - maxIndexMagnitude;
                std::get<2>(origin) = sgn(c) * (maxIndexMagnitude - 1);
                a = vxNormalized - std::get<0>(origin);
                b = vyNormalized - std::get<1>(origin);
                std::get<0>(ix) = std::get<0>(origin) + sgn(a);
                std::get<1>(iy) = std::get<1>(origin) + sgn(b);
                std::get<2>(iz) = std::get<2>(origin) - sgn(c);
                std::get<0>(e) = std::get<0>(origin) - sgn(a);
                whichExternalPoint = 0;
            }
        }
        else if(locationSum == 2) //Edge
        {
            if(location[0] == 0) //y ad z coordinate overshoot
            {
                b = vyNormalized - maxIndexMagnitude;
                std::get<1>(origin)  =  (sgn(b)) * (maxIndexMagnitude -1);
                c = vzNormalized - maxIndexMagnitude;
                std::get<2>(origin)  =  (sgn(c)) * (maxIndexMagnitude - 1);
                a = vxNormalized - std::get<0>(origin);
                std::get<0>(ix) = std::get<0>(origin) + sgn(a);
                std::get<1>(iy) = std::get<1>(origin) - sgn(b);
                std::get<2>(iz) = std::get<2>(origin) - sgn(c);
                std::get<0>(e) = std::get<0>(origin) - sgn(a);
                whichExternalPoint = 0;
            }
            else if(location[1] == 0) //z and x coordinate overshoot
            {
                c = vzNormalized - maxIndexMagnitude;
                std::get<2>(origin)  =  (sgn(c)) * (maxIndexMagnitude - 1);
                a = vxNormalized - maxIndexMagnitude;
                std::get<0>(origin)  =  (sgn(a)) * (maxIndexMagnitude - 1);
                b = vyNormalized - std::get<1>(origin);
                std::get<0>(ix) = std::get<0>(origin) - sgn(a);
                std::get<1>(iy) = std::get<1>(origin) + sgn(b);
                std::get<2>(iz) = std::get<2>(origin) - sgn(c);
                std::get<1>(e) = std::get<1>(origin) - sgn(b);
                whichExternalPoint = 1;
            }
            else //location[2] == 0 //x and y coordinate overshoot
            {
                a = vxNormalized - maxIndexMagnitude ;
                std::get<0>(origin)  =  (sgn(a)) * (maxIndexMagnitude - 1);
                b = vyNormalized - maxIndexMagnitude;
                std::get<1>(origin)  =  (sgn(b)) * (maxIndexMagnitude - 1);
                c = vyNormalized - std::get<2>(origin);
                std::get<0>(ix) = std::get<0>(origin) - sgn(a);
                std::get<1>(iy) = std::get<1>(origin) - sgn(b);
                std::get<2>(iz) = std::get<2>(origin) + sgn(c);
                std::get<2>(e) = std::get<2>(origin) - sgn(c);
                whichExternalPoint = 2;
            }
        }
        else//Corner
        {
            a = vxNormalized - maxIndexMagnitude;
            std::get<0>(origin)  =  (sgn(a)) * (maxIndexMagnitude - 1);
            b = vyNormalized - maxIndexMagnitude;
            std::get<1>(origin)  =  (sgn(b)) * (maxIndexMagnitude - 1);
            c = vzNormalized - maxIndexMagnitude;
            std::get<2>(origin)  =  (sgn(c)) * (maxIndexMagnitude - 1);
            std::get<0>(ix) = std::get<0>(origin) - sgn(a);
            std::get<1>(iy) = std::get<1>(origin) - sgn(b);
            std::get<2>(iz) = std::get<2>(origin) - sgn(c);
            whichExternalPoint = 3;
        }

        return std::array<Index, 5>{{origin, ix, iy, iz, e}};
}

std::map<Index, double> 
interpolateToGrid(const unsigned maxIndexMagnitude, const double beta, const double vx, const double vy, const double vz)
{
    const double vxNormalized = vx*beta;
    const double vyNormalized = vy*beta;
    const double vzNormalized = vz*beta;

    const auto location = pointLocation(maxIndexMagnitude, vxNormalized, vyNormalized, vzNormalized);
    unsigned locationSum = std::accumulate(location.begin(), location.end(), 0);

    std::map<Index, double> fractionalDensityChanges;
    if(locationSum == 0) //Inside vsd
    {
        const int originX =  static_cast<int>(std::nearbyint(vxNormalized)) + 5;
        const int originY =  static_cast<int>(std::nearbyint(vyNormalized)) + 5;
        const int originZ =  static_cast<int>(std::nearbyint(vzNormalized)) + 5;

        const double a = vxNormalized - originX;
        const double b = vyNormalized - originY;
        const double c = vzNormalized - originZ;

        auto const f = getFractionalDensityChanges(a, b, c, true, 0/*Dummy argument*/);
        const Index origin = std::make_tuple(originX, originY, originZ);
        const Index ix = std::make_tuple(originX + sgn(a), originY, originZ);
        const Index iy = std::make_tuple(originX, originY + sgn(b), originZ);
        const Index iz = std::make_tuple(originX, originY, originZ + sgn(c));

        const Index ex = std::make_tuple(originX - sgn(a), originY, originZ);
        const Index ey = std::make_tuple(originX, originY - sgn(b), originZ);
        const Index ez = std::make_tuple(originX, originY, originZ - sgn(c));

        const double fo = f[0];
        const double fext = f[1];
        const double fix = f[2];
        const double fiy = f[3];
        const double fiz = f[4];

        fractionalDensityChanges.emplace(origin, fo);
        fractionalDensityChanges.emplace(ix, fix);
        fractionalDensityChanges.emplace(iy, fiy);
        fractionalDensityChanges.emplace(iz, fiz);
        fractionalDensityChanges.emplace(ex, fext);
        fractionalDensityChanges.emplace(ey, fext);
        fractionalDensityChanges.emplace(ez, fext);
    }
    else
    {
        double a = 0.0;
        double b = 0.0;
        double c = 0.0;
        unsigned whichExternalPoint = 4; //0 => x, 1 => y, 2 => z, 3 => is a corner point
        auto stencil = getStencilOutsidePoint(maxIndexMagnitude, vxNormalized, vyNormalized, vzNormalized, location, a, b, c, whichExternalPoint);

        if(whichExternalPoint != 3) //Ignoring when origin is a corner point 
        {
            const Index origin = stencil[0];
            const Index ix = stencil[1];
            const Index iy = stencil[2];
            const Index iz = stencil[3];
            const Index ext = stencil[4];
            auto const f = getFractionalDensityChanges(a, b, c, false, whichExternalPoint);

            const double f0 = f[0];
            const double fext = f[1];
            const double fix = f[2];
            const double fiy = f[3];
            const double fiz = f[4];

            fractionalDensityChanges.emplace(origin, f0);
            fractionalDensityChanges.emplace(ix, fix);
            fractionalDensityChanges.emplace(iy, fiy);
            fractionalDensityChanges.emplace(iz, fiz);
            fractionalDensityChanges.emplace(ext, fext);
        }
    }
    return fractionalDensityChanges;
}

}
