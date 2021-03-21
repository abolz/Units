#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

[[nodiscard]] constexpr auto operator""_degC(long double x) noexcept {
    return DegCelsius(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_degC(unsigned long long x) noexcept {
    return DegCelsius(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_degRa(long double x) noexcept {
    return DegRankine(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_degRa(unsigned long long x) noexcept {
    return DegRankine(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_degF(long double x) noexcept {
    return DegFahrenheit(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_degF(unsigned long long x) noexcept {
    return DegFahrenheit(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_degRe(long double x) noexcept {
    return DegReaumur(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_degRe(unsigned long long x) noexcept {
    return DegReaumur(static_cast<double>(x));
}

//static void test0()
//{
//    {
//        constexpr auto t00 = Fahrenheit(Celsius(5.0)) + Fahrenheit(9.0);
//        constexpr auto t01 = Celsius(5.0) + Celsius(Fahrenheit(9.0));
//        static_assert(t00 == t01);
//    }
//}

static void test()
{
    {
        constexpr auto t00 = DegCelsius(10_K);
        constexpr auto t01 = Kelvin(t00);
        static_assert(t01.count_internal() == 10.0);
    }

    {
        constexpr auto t00 = DegCelsius(0_K);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = -273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = DegFahrenheit(0_K);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = -459.67;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(0_K);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 0.0;
        static_assert(x02 == v02);

        constexpr auto t03 = DegReaumur(0_K);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = -218.52;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = Kelvin(0_degC);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = DegFahrenheit(0_degC);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 32.0;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(0_degC);   // t_R = 9/5 t_C + 491.67
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 491.67;
        static_assert(x02 == v02); // v1 FAIL

        constexpr auto t03 = DegReaumur(0_degC);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 0.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = Kelvin(100_degC);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 100.00 + 273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = DegFahrenheit(100_degC);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 212.0;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(100_degC);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 100.0 * (9.0 / 5.0) + 491.67;
        static_assert(x02 == v02);  // v2 FAIL

        constexpr auto t03 = DegReaumur(100_degC);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 80.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = Kelvin(0_degF);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 255.37222222222222222222222222222222222222222222222222222222222222222222222222222;
        static_assert(x00 == v00);

        constexpr auto t01 = DegCelsius(0_degF);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = -17.77777777777777777777777777777777777777777777777777777777777777777777777777777;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(0_degF);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 459.67;
        static_assert(x02 == v02);

        constexpr auto t03 = DegReaumur(0_degF);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = -14.22222222222222222222222222222222222222222222222222222222222222222222222222222;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = Kelvin(32_degF);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 32.0 * (5.0 / 9.0) + (45967.0 * 5.0) / (100.0 * 9.0);
        static_assert(x00 == v00);

        constexpr auto t01 = DegCelsius(32_degF);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 32.0 * (5.0 / 9.0) - (32.0 * 5.0) / 9.0;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(32_degF);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 32.0 + 459.67;
        static_assert(x02 == v02);

        constexpr auto t03 = DegReaumur(32_degF);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 32.0 * (4.0 / 9.0) - (32.0 * 4.0) / 9.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = Kelvin(0_degRe);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 0.0 * 1.25 + 273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = DegCelsius(0_degRe);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 0.0 * 1.25;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(0_degRe);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 0.0 * 2.25 + 491.67;
        static_assert(x02 == v02);

        constexpr auto t03 = DegFahrenheit(0_degRe);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 0.0 * 2.25 + 32.0;
        static_assert(x03 == v03);
    }

    {
        constexpr auto t00 = Kelvin(80_degRe);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = 80.0 * 1.25 + 273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = DegCelsius(80_degRe);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 80.0 * 1.25;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(80_degRe);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 80.0 * 2.25 + 491.67;
        static_assert(x02 == v02);

        constexpr auto t03 = DegFahrenheit(80_degRe);
        constexpr auto x03 = t03.count_internal();
        constexpr auto v03 = 80.0 * 2.25 + 32.0;
        static_assert(x03 == v03);
    }
}

static void test2()
{
    using HeightAboveSeaLevel = Absolute<Meters>;
    using HeightAboveLocal = Absolute<Meters, Ratio<200>>;

    constexpr HeightAboveLocal h1(123);
    constexpr HeightAboveSeaLevel h2(h1);
}
