#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

using Energy
    = Joules;
using Mass
    = Kilograms;
using Velocity
    = decltype(1_m / 1_s);
using Momentum
    = decltype(1_kg * (1_m / 1_s));

static Energy total_energy(const Momentum p, const Mass m, const Velocity c)
{
    return sqrt(pow<2>(p * c) + pow<2>(m * pow<2>(c)));
}
