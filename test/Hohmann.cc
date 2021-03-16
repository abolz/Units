#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

#include <math.h>

using namespace uom;
using namespace uom::literals;

using MetersPerSecond
    = decltype(1_m / 1_s);
using SquareSeconds
    = decltype(pow<2>(1_s));
using SquareMetersPerSquareSecond
    = decltype(pow<2>(1_m) / pow<2>(1_s));
using CubicMeters
    = decltype(pow<3>(1_m));
using CubicMetersPerSquareSecond
    = decltype(pow<3>(1_m) / pow<2>(1_s));

using GravitationalParameter
    = CubicMetersPerSquareSecond;

class HohmannTransfer
{
    static MetersPerSecond dVApo(const GravitationalParameter mu, const Meters r_1, const Meters r_2)
    {
        const SquareMetersPerSquareSecond A2
            = mu / r_2;
        const Dimensionless B
            = 1.0_q - sqrt(2 / (1.0_q + r_2 / r_1));

        return sqrt(A2) * B;
    }

    static MetersPerSecond dVPer(const GravitationalParameter mu, const Meters r_1, const Meters r_2)
    {
        const SquareMetersPerSquareSecond A2
            = mu / r_1;
        const Dimensionless B
            = sqrt((2 * r_2 / r_1) / (1.0_q + r_2 / r_1)) - 1.0_q;

        return sqrt(A2) * B;
    }

public:
    const GravitationalParameter mu;
    const Meters r_1;
    const Meters r_2;
    const MetersPerSecond dV_P;
    const MetersPerSecond dV_A;

    HohmannTransfer(const GravitationalParameter _mu, const Meters _r_1, const Meters _r_2)
        : mu(_mu)
        , r_1(_r_1)
        , r_2(_r_2)
        , dV_P(dVPer(_mu, _r_1, _r_2))
        , dV_A(dVApo(_mu, _r_1, _r_2))
    {
    }
};
