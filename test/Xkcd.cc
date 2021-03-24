#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

// See:
// https://what-if.xkcd.com/11/

// 71, 73, 79
#if 1
using Birds  = Quantity<Unit<Conversion<Ratio<1>>, Kind<Dimension<71>, uom::kinds::Simple>>>;
using Poops  = Quantity<Unit<Conversion<Ratio<1>>, Kind<Dimension<73>, uom::kinds::Simple>>>;
using Mouths = Quantity<Unit<Conversion<Ratio<1>>, Kind<Dimension<79>, uom::kinds::Simple>>>;
#else
using Birds  = Quantity<Unit<Conversion<Ratio<1>>, Kind<Dimension<71>, uom::kinds::Simple>>>;
using Poops  = Quantity<Unit<Conversion<Ratio<1>>, Kind<Dimension<73>, uom::kinds::Simple>>>;
using Mouths = Quantity<Unit<Conversion<Ratio<1>>, Kind<Dimension<79>, uom::kinds::Simple>>>;
#endif

inline constexpr auto k4Pi = 4 * impl::kPi;
inline constexpr auto kEarthRadius = 6371.0_km;

inline constexpr Years answer(1 / ((Birds(300.0e+9) / (k4Pi * (kEarthRadius * kEarthRadius))) * ((Poops(1) / Birds(1)) / 1_h) * (16_h / 1_d) * (Mouths(1) / Poops(1)) * (15_cm2 / Mouths(1))));
