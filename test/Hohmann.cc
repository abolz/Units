#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

#include <math.h>

using namespace uom;
using namespace uom::literals;

using MetresPerSecond
    = decltype(1_m / 1_s);
using SquareSeconds
    = decltype(pow<2>(1_s));
using SquareMetresPerSquareSecond
    = decltype(pow<2>(1_m) / pow<2>(1_s));
using CubicMetres
    = decltype(pow<3>(1_m));
using CubicMetresPerSquareSecond
    = decltype(pow<3>(1_m) / pow<2>(1_s));

using GravitationalParameter
    = CubicMetresPerSquareSecond;

class HohmannTransfer
{
    static MetresPerSecond dVApo(const GravitationalParameter mu, const Metres r_1, const Metres r_2)
    {
        const SquareMetresPerSquareSecond A2
            = mu / r_2;
        const Dimensionless B
            = 1.0_q - sqrt(2 / (1.0_q + r_2 / r_1));

        return sqrt(A2) * B;
    }

    static MetresPerSecond dVPer(const GravitationalParameter mu, const Metres r_1, const Metres r_2)
    {
        const SquareMetresPerSquareSecond A2
            = mu / r_1;
        const Dimensionless B
            = sqrt((2 * r_2 / r_1) / (1.0_q + r_2 / r_1)) - 1.0_q;

        return sqrt(A2) * B;
    }

public:
    const GravitationalParameter mu;
    const Metres r_1;
    const Metres r_2;
    const MetresPerSecond dV_P;
    const MetresPerSecond dV_A;

    HohmannTransfer(const GravitationalParameter _mu, const Metres _r_1, const Metres _r_2)
        : mu(_mu)
        , r_1(_r_1)
        , r_2(_r_2)
        , dV_P(dVPer(_mu, _r_1, _r_2))
        , dV_A(dVApo(_mu, _r_1, _r_2))
    {
    }
};
