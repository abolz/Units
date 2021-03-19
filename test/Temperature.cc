#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

[[nodiscard]] constexpr auto operator""_degK(long double x) noexcept {
    return DegKelvin(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_degK(unsigned long long x) noexcept {
    return DegKelvin(static_cast<double>(x));
}
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

static void test()
{
    {
        constexpr auto t00 = DegCelsius(0_degK);
        constexpr auto x00 = t00.count_internal();
        constexpr auto v00 = -273.15;
        static_assert(x00 == v00);

        constexpr auto t01 = DegFahrenheit(0_degK);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = -459.67;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(0_degK);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 0.0;
        static_assert(x02 == v02);
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
    }

    {
        constexpr auto t00 = DegKelvin(0_degC);
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
    }

    {
        constexpr auto t00 = DegKelvin(100_degC);
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
    }

    {
        constexpr auto t00 = DegKelvin(0_degF);
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
    }

    {
        constexpr auto t00 = DegKelvin(32_degF);
        constexpr auto x00 = t00.count_internal();
        //constexpr auto v00 = (32.0 + 459.67) * (5.0 / 9.0);
        constexpr auto v00 = 32.0 * (5.0 / 9.0) + static_cast<double>(45967 * 5) / (100 * 9);
        static_assert(x00 == v00);

        constexpr auto t01 = DegCelsius(32_degF);
        constexpr auto x01 = t01.count_internal();
        constexpr auto v01 = 0;
        static_assert(x01 == v01);

        constexpr auto t02 = DegRankine(32_degF);
        constexpr auto x02 = t02.count_internal();
        constexpr auto v02 = 491.67;
        static_assert(x02 == v02);
    }
}