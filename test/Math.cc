#include "doctest.h"

#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

#include "FloatCompare.h"

using namespace uom;
using namespace uom::literals;

TEST_CASE("Math - sqrt")
{
    const SquareMeters m2 = 1_m2 + 1_m * 1_m;
    const Meters m = sqrt(m2);
    CHECK(count_as<Meters>(m) == 1.4142135623730951);

    const SquareCentimeters cm2 = m2;
    const Centimeters cm = sqrt(cm2);
    CHECK(count_as<Centimeters>(cm) == 141.42135623730951);

    CHECK(count_as<Meters>(cm) == 1.4142135623730951);
}

TEST_CASE("Math - trigonometric")
{
    {
        const auto x = Radians(0.0);
        const auto sin_x = sin(x);
        const auto cos_x = cos(x);
        const auto tan_x = tan(x);
        CHECK(sin_x.count_internal() == 0.0);
        CHECK(cos_x.count_internal() == 1.0);
        CHECK(tan_x.count_internal() == 0.0);
    }
    {
        const auto x = Revolutions(1.0);
        const auto sin_x = sin(Radians(x));
        const auto cos_x = cos(Radians(x));
        const auto tan_x = tan(Radians(x));
        CHECK(almost_equal(sin_x.count_internal(), 0.0));
        CHECK(almost_equal(cos_x.count_internal(), 1.0));
        CHECK(almost_equal(tan_x.count_internal(), 0.0));
    }
    {
        const auto x = Radians(1.0);
        const auto sin_x = sin(x);
        const auto cos_x = cos(x);
        const auto tan_x = tan(x);
        CHECK(sin_x.count_internal() == 0.84147098480789650);
        CHECK(cos_x.count_internal() == 0.54030230586813977);
        CHECK(tan_x.count_internal() == 1.5574077246549023);
    }
    {
        const auto x = Degrees(Radians(1.0));
        CHECK(x.count_internal() == 57.295779513082323);
#if 0
        const auto sin_x = sin(x);
        const auto cos_x = cos(x);
        const auto tan_x = tan(x);
#endif
    }
}

TEST_CASE("Math - atan2")
{
    {
        const auto x = 1_m;
        const auto y = 1_m;
        const auto phi = atan2(y, x);
        CHECK(phi == convert_to<Radians>(Degrees(45.0)));
    }
    {
        const auto x = -100_cm;
        const auto y = -1_m;
        const auto phi = atan2(y, Meters(x));
        CHECK(phi == convert_to<Radians>(Degrees(-135.0)));
    }
    {
        const auto x = 1_m;
        const auto y = -100_cm;
        const auto phi = atan2(y, Centimeters(x));
        CHECK(phi == convert_to<Radians>(Degrees(-45.0)));
    }
    //{
    //    const auto x = 1_kg;
    //    const auto y = 100_cm;
    //    const auto phi = atan2(y, x);
    //    CHECK(phi == convert_to<Radians>(Degrees(45.0)));
    //}
}

TEST_CASE("Math - round absolute")
{
    constexpr auto t00 = convert_to<DegCelsius>(0_mK);
    constexpr auto x00 = t00.count_internal();
    constexpr auto v00 = -273.15;
    static_assert(x00 == v00);

    const DegCelsius r00 = round(t00);
    CHECK(r00.count_internal() == -273.0);

    const DegCelsius r01 = round<Ratio<10>>(t00);
    CHECK(r01.count_internal() == -270.0);

    const DegCelsius r02 = round<Ratio<1,10>>(t00);
    CHECK(r02.count_internal() == -273.2);
}
