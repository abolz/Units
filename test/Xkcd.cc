#include "doctest.h"

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

TEST_CASE("Xkcd")
{
    constexpr Years answer
        = Years(
            1 / ((Gigabirds(300) / (PI4 * Rsq))
                    * ((Poops(1) / Birds(1)) / 1_h)
                    * (16_h / 1_d)
                    * (Mouths(1) / Poops(1))
                    * (15_cm2 / Mouths(1))));

    constexpr Scalar expected_answer = 193.9538756997824294;

    static_assert(answer._count_internal() == expected_answer);

    CHECK(round<Ratio<5>>(answer) == Years(195));
    CHECK(round<Weeks>(answer) == Weeks(10120));
    CHECK(round(Weeks(answer)) == Weeks(10120));
}
