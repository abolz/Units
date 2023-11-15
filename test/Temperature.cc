#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

static void test0()
{
    {
        constexpr auto t00 = 0_degC + 1_K;
        constexpr auto v00 = t00.in<DegCelsius>();
        static_assert(v00 == 1.0);
        static_assert(t00._count_internal() == 1.0);

        static_assert(1_degC - 1_K == 0_degC);
        static_assert(1_degC + 1_K == 2_degC);
        static_assert(1_K + 1_degC == 2_degC);
        static_assert(1_degC - 1_degC == 0_K);

        constexpr auto v01 = (0_degC).in<Kelvin>();
        static_assert(v01 == 273.15);

        // should not compile:
//      static_assert(1_K - 1_degC == 0_degC);
//      static_assert(1_K + 2_degC / 2 == 2_degC);
//      static_assert(1_degC + 1_degC == 2_degC);
//      static_assert(2 * 1_degC == 2_degC);
//      static_assert(1_K * 1_degC == 1_degC);

        constexpr auto v02 = (0_K).in<DegFahrenheit>();
        static_assert(v02 == -459.67);
    }
}

static void test()
{
    {
        constexpr auto t00 = (0_mK).as<DegCelsius>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = -273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = (0_K).as<DegFahrenheit>();
        constexpr auto x01 = t01._count_internal();
        constexpr auto v01 = -459.67;
        static_assert(x01 == v01);

        constexpr auto t02 = (0_K).as<DegRankine>();
        constexpr auto x02 = t02._count_internal();
        constexpr auto v02 = 0.0;
        static_assert(x02 == v02);

        constexpr auto t03 = (0_K).as<DegReaumur>();
        constexpr auto x03 = t03._count_internal();
        constexpr auto v03 = -218.52;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = (0_degC).as<Millikelvin>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = 273150.0;
        static_assert(x00 == v00);

        constexpr auto t01 = (0_degC).as<DegFahrenheit>();
        constexpr auto x01 = (0_degC).in<DegFahrenheit>();
        constexpr auto v01 = 32.0;
        static_assert(x01 == v01);

        constexpr auto t02 = (0_degC).as<DegRankine>(); // t_R = 9/5 t_C + 491.6
        constexpr auto x02 = t02._count_internal();
        constexpr auto v02 = 491.67;
        static_assert(x02 == v02);

        constexpr auto t03 = (0_degC).as<DegReaumur>();
        constexpr auto x03 = t03._count_internal();
        constexpr auto v03 = 0.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = (100_degC).as<Kelvin>();
        constexpr auto x00 = t00.in<Kelvin>();
        constexpr auto v00 = 373.15;
        static_assert(x00 == v00);

        constexpr auto t01 = (100_degC).as<DegFahrenheit>();
        constexpr auto x01 = t01._count_internal();
        constexpr auto v01 = 212.0;
        static_assert(x01 == v01);

        constexpr auto t02 = (100_degC).as<DegRankine>();
        constexpr auto x02 = t02._count_internal();
        constexpr auto v02 = 100.0 * (9.0 / 5.0) + 491.67; // ~671.67
        static_assert(x02 == v02);

        constexpr auto t03 = (100_degC).as<DegReaumur>();
        constexpr auto x03 = t03._count_internal();
        constexpr auto v03 = 80.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = (0_degF).as<Kelvin>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = 255.37222222222222222222222222222222222222222222222222222222222222222222222222222;
        static_assert(x00 == v00);

        constexpr auto t01 = (0_degF).as<DegCelsius>();
        constexpr auto x01 = t01._count_internal();
        constexpr auto v01 = -17.77777777777777777777777777777777777777777777777777777777777777777777777777777;
        static_assert(x01 == v01);

        constexpr auto t02 = (0_degF).as<DegRankine>();
        constexpr auto x02 = t02._count_internal();
        constexpr auto v02 = 459.67;
        static_assert(x02 == v02);

        constexpr auto t03 = (0_degF).as<DegReaumur>();
        constexpr auto x03 = t03._count_internal();
        constexpr auto v03 = -14.22222222222222222222222222222222222222222222222222222222222222222222222222222;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = (32_degF).as<Kelvin>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = 273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = (32_degF).as<DegCelsius>();
        constexpr auto x01 = t01._count_internal();
        constexpr auto v01 = 0.0;
        static_assert(x01 == v01);

        constexpr auto t02 = (32_degF).as<DegRankine>();
        constexpr auto x02 = t02._count_internal();
        constexpr auto v02 = 491.67;
        static_assert(x02 == v02);

        constexpr auto t03 = (32_degF).as<DegReaumur>();
        constexpr auto x03 = t03._count_internal();
        constexpr auto v03 = 0.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = (0_degRe).as<Kelvin>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = 273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = (0_degRe).as<DegCelsius>();
        constexpr auto x01 = t01._count_internal();
        constexpr auto v01 = 0.0;
        static_assert(x01 == v01);

        constexpr auto t02 = (0_degRe).as<DegRankine>();
        constexpr auto x02 = t02._count_internal();
        constexpr auto v02 = 491.67;
        static_assert(x02 == v02);

        constexpr auto t03 = (0_degRe).as<DegFahrenheit>();
        constexpr auto x03 = t03._count_internal();
        constexpr auto v03 = 32.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = (80_degRe).as<Kelvin>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = 373.15;
        static_assert(x00 == v00);

        constexpr auto t01 = (80_degRe).as<DegCelsius>();
        constexpr auto x01 = t01._count_internal();
        constexpr auto v01 = 100.0;
        static_assert(x01 == v01);

        constexpr auto t02 = (80_degRe).as<DegRankine>();
        constexpr auto x02 = t02._count_internal();
        constexpr auto v02 = 80.0 * 2.25 + 491.67; // ~671.67
        static_assert(x02 == v02);

        constexpr auto t03 = (80_degRe).as<DegFahrenheit>();
        constexpr auto x03 = t03._count_internal();
        constexpr auto v03 = 212.0;
        static_assert(x03 == v03);
    }
}

template <typename R, typename A>
using ScaledAbsolute
    = typename Absolute<
        typename ScaledQuantity<Conversion<R>, typename A::relative_type>::type,
        typename std::ratio_divide<typename A::zero, R>::type>::type;

static void test2()
{
    using MilliDegCelsius = ScaledAbsolute<Ratio<1, 1000>, DegCelsius>;
    // using MilliDegCelsius = Absolute<Millikelvin, Ratio<27315000, 100>>;

    {
        constexpr auto t00 = (0_mK).as<DegCelsius>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = -273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = (t00).as<MilliDegCelsius>();
        constexpr auto x01 = t01._count_internal();
        constexpr auto v01 = -273150.0;
        static_assert(x01 == v01);
    }
    {
        constexpr auto t00 = (MilliDegCelsius(1000.0)).as<DegCelsius>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = 1.0;
        static_assert(x00 == v00);
    }
    {
        constexpr auto t00 = (1_degC).as<MilliDegCelsius>();
        constexpr auto x00 = t00._count_internal();
        constexpr auto v00 = 1000.0;
        static_assert(x00 == v00);
    }
}
