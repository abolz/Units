#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

static void test0()
{
    {
        constexpr auto t00 = 0_degC + 1_K;
        constexpr auto v00 = t00.count<DegCelsius>();
        static_assert(v00 == 1.0);
        static_assert(t00.count_internal() == 1.0);

        static_assert(1_degC - 1_K == 0_degC);
        static_assert(1_degC + 1_K == 2_degC);
        static_assert(1_K + 1_degC == 2_degC);
        static_assert(1_degC - 1_degC == 0_K);

        constexpr auto v01 = (0_degC).count<Kelvin>();
        static_assert(v01 == 273.15);

        // should not compile:
//      static_assert(1_K - 1_degC == 0_degC);
//      static_assert(1_K + 2_degC / 2 == 2_degC);
//      static_assert(1_degC + 1_degC == 2_degC);
//      static_assert(2 * 1_degC == 2_degC);
//      static_assert(1_K * 1_degC == 1_degC);

        constexpr auto v02 = (0_K).count<DegFahrenheit>();
        static_assert(v02 == -459.67);
    }
}

static void test()
{
    {
        constexpr auto t00 = convert_to<DegCelsius>(0_mK);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = -273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = convert_to<DegFahrenheit>(0_K);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = -459.67;
        static_assert(x01 == v01);

        constexpr auto t02 = convert_to<DegRankine>(0_K);
        constexpr auto x02 = t02.count<Rankine>();
        constexpr auto v02 = 0.0;
        static_assert(x02 == v02);

        constexpr auto t03 = convert_to<DegReaumur>(0_K);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = -218.52;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = convert_to<Millikelvin>(0_degC);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 273150.0;
        static_assert(x00 == v00);

        constexpr auto t01 = convert_to<DegFahrenheit>(0_degC);
        constexpr auto x01 = (0_degC).count<DegFahrenheit>(); // t01.count_internal();
        constexpr auto v01 = 32.0;
        static_assert(x01 == v01);

        constexpr auto t02 = convert_to<DegRankine>(0_degC);   // t_R = 9/5 t_C + 491.67
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 491.67;
        static_assert(x02 == v02);

        constexpr auto t03 = convert_to<DegReaumur>(0_degC);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 0.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = convert_to<Kelvin>(100_degC);
        constexpr auto x00 = t00.count<Kelvin>();
        constexpr auto v00 = 373.15;
        static_assert(x00 == v00);

        constexpr auto t01 = convert_to<DegFahrenheit>(100_degC);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 212.0;
        static_assert(x01 == v01);

        constexpr auto t02 = convert_to<DegRankine>(100_degC);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 100.0 * (9.0 / 5.0) + 491.67; // ~671.67
        static_assert(x02 == v02);

        constexpr auto t03 = convert_to<DegReaumur>(100_degC);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 80.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = convert_to<Kelvin>(0_degF);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 255.37222222222222222222222222222222222222222222222222222222222222222222222222222;
        static_assert(x00 == v00);

        constexpr auto t01 = convert_to<DegCelsius>(0_degF);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = -17.77777777777777777777777777777777777777777777777777777777777777777777777777777;
        static_assert(x01 == v01);

        constexpr auto t02 = convert_to<DegRankine>(0_degF);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 459.67;
        static_assert(x02 == v02);

        constexpr auto t03 = convert_to<DegReaumur>(0_degF);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = -14.22222222222222222222222222222222222222222222222222222222222222222222222222222;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = convert_to<Kelvin>(32_degF);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = convert_to<DegCelsius>(32_degF);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 0.0;
        static_assert(x01 == v01);

        constexpr auto t02 = convert_to<DegRankine>(32_degF);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 491.67;
        static_assert(x02 == v02);

        constexpr auto t03 = convert_to<DegReaumur>(32_degF);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 0.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = convert_to<Kelvin>(0_degRe);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = convert_to<DegCelsius>(0_degRe);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 0.0;
        static_assert(x01 == v01);

        constexpr auto t02 = convert_to<DegRankine>(0_degRe);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 491.67;
        static_assert(x02 == v02);

        constexpr auto t03 = convert_to<DegFahrenheit>(0_degRe);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 32.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = convert_to<Kelvin>(80_degRe);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 373.15;
        static_assert(x00 == v00);

        constexpr auto t01 = convert_to<DegCelsius>(80_degRe);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 100.0;
        static_assert(x01 == v01);

        constexpr auto t02 = convert_to<DegRankine>(80_degRe);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 80.0 * 2.25 + 491.67; // ~671.67
        static_assert(x02 == v02);

        constexpr auto t03 = convert_to<DegFahrenheit>(80_degRe);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 212.0;
        static_assert(x03 == v03);
    }
}

static void test2()
{
    using GlobalHeight = Meters; // Absolute<Meters>;
    using LocalHeight = Absolute<Meters, Ratio<200>>;

    constexpr LocalHeight h1(123.0);
    static_assert(h1.count_internal() == 123.0);
    static_assert(h1.count<GlobalHeight>() == 323.0);

    constexpr GlobalHeight h2 = convert_to<GlobalHeight>(h1);
    static_assert(h2.count_internal() == 323.0);
    static_assert(h2.count<LocalHeight>() == 123.0);
}
