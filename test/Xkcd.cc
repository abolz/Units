#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

// See:
// https://what-if.xkcd.com/11/

using Birds  = TaggedQuantity<Entities, struct _birds>;
using Poops  = TaggedQuantity<Entities, struct _poops>;
using Mouths = TaggedQuantity<Entities, struct _mouths>;

using Gigabirds = ScaledQuantity<Conversion<Ratio<1'000'000'000>>, Birds>;

static constexpr auto PI4 = 4 * 3.1415;
static constexpr auto R = 6371_km;
static constexpr auto Rsq = R * R;

static constexpr Years answer
    = Years(
        1 / ((Gigabirds(300.0) / (PI4 * Rsq))
                * ((Poops(1_ent) / Birds(1_ent)) / 1_h)
                * (16_h / 1_d)
                * (Mouths(1_ent) / Poops(1_ent))
                * (15_cm2 / Mouths(1_ent))));

static_assert(answer.count_internal() == 193.9538756997824294);
//static_assert(constexpr_math::ieee754_round(answer.count_internal()) == 194);
//static_assert(round(answer).count_internal() == 194);

//int main()
//{
//    const auto y1 = round<Ratio<5>>(answer); // round to 5-years
//    const auto y2 = round<Weeks>(answer);
//    const auto y3 = round(Weeks(answer));
//    return 0;
//}
