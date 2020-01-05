#if 0
#include "benri/si/si.h"
#include "benri/cmath.h"
#include "benri/quantity.h"
#include "benri/quantity_point.h"

using namespace benri;
using namespace benri::si;    //Import si literals.
using namespace benri::casts; //Import casts into namespace for ADL to work.

//static void test001()
//{
//    constexpr auto l1 = 1_centi * metre;
//    constexpr auto l2 = 1_metre;
//    constexpr auto l3 = l1 + l2;
//}

static constexpr void test002(quantity<metre_t> m)
{
}

static void test003()
{
    test002(1_metre);
    test002(1_kilo * metre);
}

int main()
{
    //Set size of the cake.
    constexpr auto height = 4_centi * metre;
    constexpr auto diameter = 30_centi * metre;

    //Calculate the volume.
    constexpr auto volume = height * constant::pi * square(diameter / 2.0);

    //Set the density of the batter.
    constexpr auto density = 1_gram / cubic(centi * metre);

    //Calculate the mass of the cake.
    constexpr auto mass = density * volume;

    ////Print the recipe.
    //std::cout
    //    << "  vanilla cake  \n"
    //    << "----------------\n"
    //    << ceil(mass * 0.028).value()  << " gram butter\n"
    //    << ceil(mass * 0.099).value()  << " gram sugar\n"
    //    << ceil(mass * 0.0635).value() << " gram flour\n"
    //    << ceil(mass / 2828.0).value() << " pack of backing soda\n"
    //    << floor(mass * 0.002).value() << " eggs\n"
    //    << ceil(mass * 0.1765).value()  << " gram quark cream\n"
    //    << ceil(mass * 0.0705).value() << " gram oil\n"
    //    << ceil(mass / 2828.0).value() << " pack of vanille aroma\n"
    //    << ceil(mass / 2828.0).value() << "/2 litre milk\n"
    //    << std::flush;
}
#endif

#if 0
#include "../src/FsUnits.h"

#include <limits>
#include <chrono>

using namespace DkCore;
using namespace DkCore::Literals;

template <typename ...> struct Incomplet;

static constexpr double Pi = DkCore::Impl::PowerOfPi<1>::value;

int main()
{
    constexpr std::chrono::seconds chr1 = std::chrono::nanoseconds{1};

    //constexpr auto m0 = 1 * Meter{};
    constexpr Metres<> m1 = Millimetres<>{1};
    constexpr Millimetres<> m2 = Metres<>{1};
    constexpr auto m3 = m1 + m2;
    constexpr auto m4 = 1_m + 1_cm + 1_mm;
    //constexpr auto m5 = 1 / (1 / 2.98873687628376_m);
    //static_assert(m5 == 2.98873687628376_m, "");

    constexpr auto m6 = 1_mi + 1_yd + 1_ft + 1_in;
    constexpr auto m7 = Micrometres{m6};
    constexpr auto m8 = 1_mi + 1_yd;
    constexpr auto m9 = 1_yd + 1_ft;
    constexpr auto m10 = 1_ft + 1_in;
    constexpr auto m11 = 1_in + 1_cm;
    constexpr auto cm11 = Centimetres{m11};
    constexpr auto cm11x = 1.0 + 2.54;
    constexpr auto m12 = Centimetres{1_in};
    constexpr Centimetres m13 = 1_in;
    constexpr auto m14 = 1_in;
    //using InchType = units::Inch::conversion;
    //Incomplet<InchType>{};
    //Incomplet<decltype(m11)::conversion>{};
    constexpr auto m15 = Centimetres{Metres{Yards{Inches{1_cm}}}};
    //Incomplet<decltype(m15)::conversion>{};
    constexpr auto m16 = Centimetres{Inches{1_cm}};
    //constexpr auto m17 = Centimetres{Yards{Inches{1_cm}}};
    constexpr auto m18 = 3.2347687234_mi * 1_cm * 2.2347326487236487234_in / 1e+2_yd * 1e+2_yd / 2.2347326487236487234_in / 3.2347687234_mi;
    constexpr auto m18cm = Centimetres{m18};
    static_assert(m18cm == 1_cm, "");
    constexpr auto m19 = 3.2347687234_mi + 1_cm + 2.2347326487236487234_in - 2.2347326487236487234_in - 3.2347687234_mi + 1_cm;
    //Incomplet<decltype(m19)::conversion>{};
    static_assert(m19 == 2_cm, "");

    constexpr Radians r1 = 90_deg;
    constexpr auto r2 = Degrees{Radians{90_deg}};
    //Incomplet<decltype(r2)::conversion>{};
    constexpr auto r22 = Radians{Degrees{1_rad}};
    //Incomplet<decltype(r22)::conversion>{};
    constexpr auto r3 = 1_rad + 1_deg;
    constexpr auto r4 = 1_deg + 1_rad;
    constexpr auto r5 = 1_deg + 2_deg / 2 + 3_gon / 3;
    constexpr auto r6 = 1_gon + 1_deg;
    constexpr auto r7 = 1_deg + 1_gon;
    // constexpr Radians r8 = 1_rev;
    // constexpr Degrees r9 = r8;
    // static_assert(r8 == Radians(2 * Pi), "");

    constexpr auto v1 = 1_m / 1_s;
    constexpr auto v2 = 3_km / 1_h;
    constexpr auto v3 = v1 + v2;

    constexpr auto N1 = 1_kg * 1_m / (1_s * 1_s);
    constexpr auto N2 = 1_m / 1_s / 1_s * 1_kg;
    constexpr auto N3 = 1_kg / (1_s * 1_s) * 1_m;
    static_assert(std::is_same<decltype(N1)::conversion, decltype(N2)::conversion>::value, "");
    static_assert(std::is_same<decltype(N1)::conversion, decltype(N3)::conversion>::value, "");
    //Incomplet<decltype(N1)>{};
    //Incomplet<decltype(N1)::conversion>{};
    //Incomplet<decltype(N2)::conversion>{};
    constexpr auto N4 = 1_g * 1_cm / (1_ms * 1_ms);
    //Incomplet<decltype(N4)>{};

    //constexpr auto xxx1 = 1_m + 1_s;

    //constexpr auto ouch1 = 1_m / 1_m + 1_s / 1_s + ((1_m/1_s) / (1_m/1_s)) + Pi;
    //static_assert(ouch1 == 3 + Pi, "");

    constexpr auto test1 = 1_m;
    constexpr auto test2 = 1_s;
    constexpr auto test3 = 1_m / 1_s + (1 / 1_s) * 1_m;
    constexpr auto test4 = 1_m / 1_s + 1_km / 1_h;
}
#endif

#if 1
#include "../src/Units.h"

#include <stdio.h>
#include <chrono>

template <typename...> struct Incomplet;

//namespace sc {
//    namespace kinds {
//        struct Rainfall : UnitKind<Rainfall, decltype(Volume::dimension{} / Area::dimension{})> {};
//
//        //constexpr auto operator/(Volume, Area) { return Rainfall{}; }
//    }
//
//    using Rainfall  = Unit< Ratio<1, 1>, kinds::Rainfall >;
//    using Rainfalls = Quantity< Rainfall >;
//}
//
//namespace sc {
//    namespace kinds {
//        struct LengthTime : UnitKind<LengthTime, decltype(Length::dimension{} * Time::dimension{})> {};
//
//        constexpr auto operator*(Time, Length) { return LengthTime{}; }
//        constexpr auto operator*(Length, Time) { return LengthTime{}; }
//        constexpr auto operator/(LengthTime, Time) { return Length{}; }
//        constexpr auto operator/(LengthTime, Length) { return Time{}; }
//        constexpr auto operator/(LengthTime, LengthTime) { return Dimensionless{}; }
//    }
//
//    using MeterSecond = Unit< Ratio< 1, 1>, kinds::LengthTime>;
//    using MeterSeconds = Quantity<MeterSecond>;
//}

#if RATIO_ROOT()
namespace sc {
    using impl::Root;
    using impl::Power;
    static_assert(Root(1, 1) == 1, "");
    static_assert(Root(1, 2) == 1, "");
    static_assert(Root(1, 10) == 1, "");
    static_assert(Root(1, INTMAX_MAX) == 1, "");
    static_assert(Root(1, INTMAX_MAX / 2) == 1, "");
    static_assert(Root(2, INTMAX_MAX) == 1, "");
    static_assert(Root(2, INTMAX_MAX / 2) == 1, "");
    static_assert(Root(2, 2) == 1, "");
    static_assert(Root(2, 1) == 2, "");
    static_assert(Root(3, 2) == 1, "");
    static_assert(Root(4, 2) == 2, "");
    static_assert(Root(5, 2) == 2, "");
    static_assert(Root(8, 2) == 2, "");
    static_assert(Root(9, 2) == 3, "");
    static_assert(Root(10, 2) == 3, "");
    static_assert(Root(15, 2) == 3, "");
    static_assert(Root(16, 2) == 4, "");
    static_assert(Root(27, 3) == 3, "");
    static_assert(Root(152399025, 2) == 12345, "");
    static_assert(Root(1881365963625, 3) == 12345, "");
    static_assert(Root(23225462820950625, 4) == 12345, "");
    static_assert(Root(3530945043777457216, 6) == 1234, "");
    static_assert(Root(8650415919381337933, 17) == 13, "");
    static_assert(Root(8650415919381337934, 17) == 13, "");
    static_assert(Root(INTMAX_MAX, 17) == 13, "");
    static_assert(Root(INTMAX_MAX, INTMAX_MAX) == 1, "");
    static_assert(Root(INTMAX_MAX, 2) == 3037000499, "");
    static_assert(Root(INTMAX_MAX, 3) == 2097151, "");
    static_assert(Root(INTMAX_MAX, 4) == 55108, "");
    static_assert(Root(INTMAX_MAX, 1) == INTMAX_MAX);
}
#endif

//template <typename C1, typename C2>
//constexpr auto avg_speed(sc::Length<C1> l, sc::Time<C2> t) noexcept
//{
//    return l / t;
//    //return sc::Velocity{l / t};
//}

using namespace sc;
using namespace sc::literals;

#if 1
template <typename C>
using Velocity = Quantity< Unit< C, kinds::Velocity > >;
#else
template <typename C>
using Velocity = Quantity< Unit< C, DivKinds<kinds::Length, kinds::Time> > >;
#endif

template <typename C>
static void takesVelocity(Velocity<C> /*v*/)
{
}

static void test001()
{
    //constexpr auto x1 = 1_mps;
    //constexpr auto y1 = kind_cast<kinds::Velocity>(1_m / 1_s);
    ////constexpr auto y1 = 1_m / 1_s;
    //static_assert(x1 == y1);

    constexpr auto x2 = 1_mps;
    constexpr auto y2 = 1_m / 1_s;
    //static_assert(x2 == y2);
    takesVelocity(x2);
    takesVelocity(y2);
}

template <typename L, typename R> using Add = decltype(std::declval<L>() + std::declval<R>());
template <typename L, typename R> using Sub = decltype(std::declval<L>() - std::declval<R>());
template <typename L, typename R> using Mul = decltype(std::declval<L>() * std::declval<R>());
template <typename L, typename R> using Div = decltype(std::declval<L>() / std::declval<R>());

template <typename From, typename To>
using ConvertFromTo = decltype(std::declval<From>().convert_to(std::declval<To>()));

namespace details
{
    template <typename T>
    struct Void
    {
        using type = void;
    };

    template <typename T>
    using Void_t = typename Void<T>::type;

    template <typename AlwaysVoid, template <typename...> class Op, typename... Args>
    struct Detector
    {
        static constexpr bool value = false;
    };

    template <template <typename...> class Op, typename... Args>
    struct Detector<Void_t<Op<Args...>>, Op, Args...>
    {
        static constexpr bool value = true;
    };
}

template <template <typename...> class Op, typename... Args>
inline constexpr bool Compiles = details::Detector<void, Op, Args...>::value;

static void test002()
{
    {
//#if !UNITS_IGNORE_KIND()
//        constexpr auto x = 1_m / 1_m;
//        constexpr auto y = 1_s / 1_s;
//        constexpr auto z = x * y;
//        static_assert(!Compiles<Add, decltype(x), decltype(y)>, "");
//#endif
    }
    {
        constexpr auto x = 1_m;
        constexpr auto y = 1_s;
        //constexpr auto z = x + y;
        static_assert( Compiles<Mul, decltype(x), decltype(y)>, "");
        static_assert(!Compiles<Add, decltype(x), decltype(y)>, "");
    }
    {
        constexpr auto m01 = 1_m * 1_s;
        constexpr auto m02 = 1_s * 1_m;
        //static_assert(std::is_same_v< decltype(m01)::kind, kinds::Product<kinds::Length, kinds::Time> >, "");
        //static_assert(kinds::CompatibleKinds<decltype(m01)::kind, decltype(m02)::kind>::value, "");
        //constexpr auto m03 = m01 + m02;
    }
    //{
    //    constexpr auto t01 = 1_m / 1_s;
    //    constexpr auto t02 = 1_s / 1_m;
    //    constexpr auto t03 = t01 + t02;
    //}
}

static void test003()
{
//#if UNITS_IGNORE_KIND()
//    constexpr auto t0 = 1_rad;
//    constexpr auto t1 = 1_sr;
//    constexpr auto t2 = t0 + Radians{t1};
//    //constexpr auto t3 = Radians{1_m};
//    //constexpr auto t3 = t0 + t1;
//#endif
}

static void test004()
{
    constexpr auto t0 = 1_m / 1_s;
    constexpr auto t1 = 1_km / 1_h;
    constexpr auto t2 = 1_mps;
    constexpr auto t3 = 1_kmph;
    static_assert(std::is_convertible_v<decltype(t0), MetresPerSecond>, "");
    static_assert(std::is_convertible_v<decltype(t2), MetresPerSecond>, "");
}

//#if UNITS_HAS_ANY()
//static constexpr int takesLength(Centimetres)
//{
//    return 0;
//}
//
//static void testFlatten()
//{
//    constexpr auto len1 = 1_m;
//    constexpr auto dur1 = 1_s;
//    //constexpr auto xxx1 = len1 * dur1 / dur1;
//    constexpr auto xxx1 = len1 / dur1 * dur1;
//    //constexpr auto res1 = takesLength(xxx1);
//    constexpr auto res2 = takesLength(flatten(xxx1));
//    //constexpr auto res3 = takesLength(as_any(xxx1));
//    //constexpr auto res4 = takesLength(as_length(xxx1));
//    //constexpr auto res5 = takesLength(as_length(1_s));
//    //constexpr auto xxx2 = kind_cast<kinds::Length>(xxx1);
//    //constexpr auto res6 = takesLength(xxx2);
//    //constexpr auto res6 = takesLength(as_any_quantity(1_m)); // any_kind_cast(1_s));
//    constexpr auto res7 = takesLength(1_m);
//}
//#endif

template <typename C>
using AngularVelocity = Quantity< Unit< C, kinds::AngularVelocity > >;

template <typename C>
static void takesAngularVelocity(AngularVelocity<C>)
{
}

static void test009()
{
    takesAngularVelocity(1_rad / 1_h);
    //takesAngularVelocity(1 / 1_s);
}

#if UNITS_HAS_MATH()
static void testFma0()
{
    constexpr auto sum0 = 3_mm * 1_s + 2_km * 1_s;
    //Incomplet<decltype(sum0)>{};
    constexpr auto sum1 = Fma(3_mm, 1_s, 2_km * 1_s);
    //Incomplet<decltype(sum1)>{};
    static_assert(sum0.count() == sum1.count());
}

static void testFma1()
{
    constexpr auto sum0 = 1_mm * 1_h + 2_km * 1_s;
    //Incomplet<decltype(sum0)>{};
    //constexpr auto sum1 = Fma(1_mm, 1_h, 2_km * 1_s);
    //Incomplet<decltype(sum1)>{};
    //static_assert(sum0.count() == sum1.count());
}

static void testFma2a()
{
    constexpr auto sum0 = 1_m * (1 / 1_s) + 1_m * (1 / 1_s);
    //Incomplet<decltype(sum0)>{};
    constexpr auto sum1 = Fma(1_m, 1 / 1_s, 1_m * (1 / 1_s));
    //Incomplet<decltype(sum1)>{};
    static_assert(sum0.count() == sum1.count());
}

static void testFma2()
{
    constexpr auto sum0 = 1_km * (1 / 1_h) + 1_m * (1 / 1_s);
    //Incomplet<decltype(sum0)>{};
    //constexpr auto sum1 = Fma(1_km, 1 / 1_h, 1_m * (1 / 1_s));
    //Incomplet<decltype(sum1)>{};
    //static_assert(sum0.count() == sum1.count());
}
#endif

#if 1
static void test1()
{
    //constexpr Metres zzz = 1_s;
    constexpr Centimetres zzz = 1_m;
    constexpr Centimetres zzz2{1_m};
    //constexpr Centimetres zzz3{1_s};
    //constexpr Metres z2 = Seconds{1}.convert_to(Metres{});

    ////auto vel1 = 1_mps;
    ////auto vel1a = 1_m * (1 / 1_s);
    ////auto vel1b = vel1 + kind_cast<kinds::Velocity>(vel1a);
    ////auto vel2 = 1_kmph;
    //////vel1 = vel2;
    ////vel1 = MetresPerSecond{vel2};

    //auto len1 = 1_m;
    //len1 = 1_mm;
    auto len2 = 1_mm;
    len2 = 1_m;

    //constexpr Quantity m = 1_m;
    //constexpr Quantity s = 1_s;
    //constexpr Quantity z = m + s;

    auto vel01 = 1_kmph;
    //KilometresPerHour vel02 = 1_mps;
    KilometresPerHour vel03 = KilometresPerHour{1_mps};

    constexpr auto phi0 = 1_gon + 1_deg;
    //Incomplet<decltype(phi0)>{};
    constexpr auto phi1 = Degrees{phi0};
    //constexpr auto phi2 = 1 / 1_gon + 1 / 1_deg;
    constexpr auto phi3 = 1 / 1_gon;
    //constexpr auto val3 = phi3.value();
    constexpr auto phi4 = 1 / 1_deg;
    //constexpr auto val4 = phi4.value();
    constexpr auto phi5 = phi3 + phi4;
    //constexpr auto val5 = phi5.value();
    //Incomplet<decltype(phi5)>{};

    static_assert(!Compiles<Add, Radians, Degrees>, "");
    //Incomplet<decltype(phi6)>{};
    constexpr auto phi6a = Radians{1_deg} + 1_rad;
    //Incomplet<decltype(phi6a)>{};
    constexpr auto phi6b = 1_deg + Degrees{1_rad};
    //Incomplet<decltype(phi6b)>{};

    constexpr auto phi7_0 = 1_deg + 1_rev;
    //Incomplet<decltype(phi7_0)>{};
    constexpr auto phi8_0 = 1_deg * 1_rev;
    //Incomplet<decltype(phi8_0)>{};
    constexpr auto phi9_0 = phi7_0 / phi8_0;
    //Incomplet<decltype(phi9_0)>{};
    constexpr auto phi10_0 = 1 / 1_deg + 1 / 1_rev;
    //Incomplet<decltype(phi10_0)>{};

    constexpr auto phi7 = 1_m + 1_km;
    constexpr auto phi8 = 1_m * 1_km;
    constexpr auto phi9 = phi7 / phi8;
    //Incomplet<decltype(phi9)>{};
    constexpr auto phi10 = 1 / 1_m + 1 / 1_km;
//#if !UNITS_IGNORE_KIND()
//    static_assert(!Compiles<Add, decltype(phi9), decltype(phi10)>, "");
//#endif
    using Phi9 = std::remove_const_t<decltype(phi9)>;
    constexpr auto phi12 = phi9 - Phi9(phi10);
    //Incomplet<decltype(phi12)>{};
//#if UNITS_IGNORE_KIND()
//    constexpr auto phi12a = phi9 - phi10;
//#endif
//#if UNITS_HAS_ANY()
//    constexpr auto phi13 = flatten(phi9) - flatten(phi10);
//    //Incomplet<decltype(phi13)>{};
//
//    constexpr Phi9 phi14 = flatten(phi10);
//
//    //constexpr auto phi15 = flatten(phi9) - phi10;
//    //Incomplet<decltype(phi15)>{};
//
//    constexpr auto phi16 = flatten(flatten(phi9) * flatten(phi10));
//    //Incomplet<decltype(phi16)>{};
//#endif
}
#endif

static void test999()
{
    constexpr auto bits01 = Bits{1};
    constexpr auto bits02 = Bits{1} + Bytes{1} + Kilobytes{1};
    //constexpr auto bits03 = Bits{1_rad};

    constexpr auto bits04 = 1_MB / 1_s;
    constexpr auto bits05 = bits04.convert_to(Gigabytes{1} / Hours{1});
    constexpr auto bits15 = bits04.convert_to(1_GB / 1_h);
    constexpr auto bits06 = 1_GB / 1_h;
    constexpr auto bits07 = bits06.convert_to(Megabytes{1} / Seconds{1});
    constexpr auto bits08 = bits06.convert_to(1_MB / 1_s);
    constexpr auto bits09 = bits06.convert_to(units::Megabyte{} / units::Second{});
}

static void test998()
{
    constexpr auto t0 = 1_ms;
    constexpr auto t1 = remove_conversion(t0);
}

static constexpr void test997(Centimetres cm)
{
}

template <typename T>
using CallTest997 = decltype(test997(std::declval<T>()));

static void test996()
{
    test997(1_cm);
    test997(1_m);
    //test997(Millimetres(1).convert_to(1_m));
    //test997(1_mm.convert_to(1_m));
    static_assert(!Compiles<CallTest997, decltype(1_mm)>, "");
}

// f(x) = x^2
static constexpr SquareMetres square(Metres x) {
    return x * x - SquareMetres(2);
}

// f'(x) = 2x
static constexpr Metres d_square(Metres x) {
    //return 2 * x;

    constexpr auto h = Metres(std::numeric_limits<double>::epsilon() * 1000);
    return Metres{ (square(x + h) - square(x - h)) / (2 * h) };
}

// x_{n+1} = x_n - f(x_n) / f'(x_n)
static constexpr Metres newton_step(Metres x) {
    return x - Metres{square(x) / d_square(x)};
}

static void testNewton()
{
    constexpr auto s1 = newton_step(2_m);
    constexpr auto s2 = newton_step(s1);
    constexpr auto s3 = newton_step(s2);
    constexpr auto s4 = newton_step(s3);
    constexpr auto s5 = newton_step(s4);
    constexpr auto s6 = newton_step(s5);
    constexpr auto s7 = newton_step(s6);
    constexpr auto s8 = newton_step(s7);
    constexpr auto s9 = newton_step(s8);
}

static constexpr Metres integrate_step(Metres x, Metres h, int n)
{
    return n <= 0
        ? 0_m
        : d_square(x) + integrate_step(x + h, h, n - 1);
}

static constexpr auto integrate(Metres a, Metres b, int n)
{
    const auto h = (b - a) / n;
    return (0.5 * h) * (d_square(a) + 2 * integrate_step(a + h, h, n - 1) + d_square(b));
}

static void testIntegrate()
{
    constexpr auto i1 = integrate(0_m, 2_m, 1);
    constexpr auto i2 = integrate(0_m, 2_m, 2);
    constexpr auto i3 = integrate(1_m, 2_m, 3);
    constexpr auto i4 = integrate(1_m, 2_m, 4);
    constexpr auto i5 = integrate(2_m, 3_m, 5);
    constexpr auto i6 = integrate(2_m, 3_m, 6);
}

static void test800()
{
    // C = K + 273.15
    // K = C - 273.15

    // F -> C:  Kelvin{x - FahrenheitZero} + CelsiusZero
    // C -> F:  Kelvin{x - CelsiusZero} + FahrenheitZero

    ////constexpr auto c01 = 20_degC;
    ////constexpr auto c02 = c01 + 100_K;
    ////constexpr auto c03 = 200_degC - 10_degC;
    ////constexpr auto c04 = Celsius{200_degC - 10_degC};
    //////constexpr auto c05 = 200_degC + 10_degC;
    ////constexpr auto c06 = 0_degC - CelsiusZero; // K
    ////constexpr auto c07 = 0_K + CelsiusZero; // Torsor<K> = C
    ////constexpr auto c08 = Kelvin{10_degF - FahrenheitZero};
    ////constexpr auto c09 = Celsius{c08 + CelsiusZero};
    ////constexpr auto c10 = Fahrenheit{ (Celsius{c09} - CelsiusZero) }; // XXXXXXXXXXXXXXXXXXXXXXXXXX this is wrong...
    ////constexpr auto c11 = Fahrenheit{ FahrenheitZero + (Celsius{c09} - CelsiusZero) };
}

//using HorMPS = Quantity< Unit< Conversion_t<1>, TaggedKind< struct HorVel_t, kinds::Velocity > > >;
//using VerMPS = Quantity< Unit< Conversion_t<1>, TaggedKind< struct VerVel_t, kinds::Velocity > > >;

//using HorMPS = Quantity< TaggedUnit< struct HorVel_t, units::MetrePerSecond > >;
//using VerMPS = Quantity< TaggedUnit< struct VerVel_t, units::MetrePerSecond > >;

using HorMPS = TaggedQuantity< struct HorVel, MetresPerSecond >;
using VerMPS = TaggedQuantity< struct VerVel, MetresPerSecond >;

static void test700_fun(HorMPS v1, VerMPS v2)
{
    //const auto v9 = v1 + v2;
    const auto v0 = v1 + HorMPS{v2};
    const auto v3 = v1 * v2;
    const auto v4 = v2 * v1;
    //const auto v5 = v3 + v4;
}

static void test700()
{
    constexpr HorMPS v1{1.0};
    constexpr VerMPS v2{1.0};
    constexpr HorMPS v3{1_mps};
    constexpr HorMPS v4{1_kmph};
    //v1 = v2;
    //v1 = HorVel{v2};
    static_assert(!Compiles<Add, HorMPS, VerMPS>, "");
    static_assert(!Compiles<Add, VerMPS, HorMPS>, "");
    //const auto v3 = v1 + v2;
    //test700_fun(v2, v1);
    test700_fun(v1, v2);
    test700_fun(HorMPS{1_m / 1_s}, VerMPS{1_m / 1_s});
}

static void test699()
{
#if UNITS_DIMENSIONLESS_ARITHMETIC()
    constexpr auto r1 = 1_m / 1_m;
    constexpr auto r2 = 1_s / 1_s;
    static_assert(Compiles<Add, decltype(r1), double>, "");
    static_assert(Compiles<Sub, decltype(r1), double>, "");
    static_assert(Compiles<Add, double, decltype(r1)>, "");
    static_assert(Compiles<Sub, double, decltype(r1)>, "");
#else
    constexpr auto r1 = 1_m / 1_m;
    constexpr auto r2 = 1_s / 1_s;
    static_assert(!Compiles<Add, decltype(r1), double>, "");
    static_assert(!Compiles<Sub, decltype(r1), double>, "");
    static_assert(!Compiles<Add, double, decltype(r1)>, "");
    static_assert(!Compiles<Sub, double, decltype(r1)>, "");
#endif
}

//struct Width : Metres {};
//struct Height : Metres {};
using Width = TaggedQuantity< struct Width_t, Metres >;
using Height = TaggedQuantity< struct Height_t, Metres >;

static constexpr void test698a(Width w, Height h)
{
    //const auto z = 1_m + 1_s;
    //const auto x = w + h;
    const auto y = w * h;
}

static void test698()
{
}

int main()
{
    constexpr auto hhhh = 1_h + 1_s;
    constexpr Metres m = 1_m;
    constexpr Metres m2222 = Metres{1_mm};
    constexpr Millimetres mm = 1_mm;
    constexpr Millimetres mmmm22222 = 1_m;
    constexpr auto mmsq = mm * mm;
    constexpr auto mm01 = mm.count();
    //constexpr auto mm02 = mm.value();
    constexpr auto m2 = 1_m + 1_cm + 1_mm + 1_nm;
    constexpr auto one_m = 1_m;
    constexpr auto two_yd = 2_yd;
    constexpr auto xxx = 1_m + 2_yd + 1_mm;
    constexpr auto xxx2 = 2_yd + 1_m + 1_mm;
    constexpr Inches in = 1_mi;
    //constexpr Miles mi = 1_in;
    static_assert(xxx == xxx2);
    using TTT = decltype(xxx2);
    //constexpr auto ttt = TTT::conv;
    //constexpr auto tttNum = TTT::conv_type::num;
    //constexpr auto tttDen = TTT::conv_type::den;
    //constexpr double xxx_value = xxx.value();
    constexpr auto yyy = quantity_cast<Millimetres>(xxx);
    constexpr auto yyy1 = Millimetres{xxx};
    //constexpr Millimetres yyy2 = xxx;
    constexpr Millimetres yyy3 = 1_m;
    static_assert(std::is_same<decltype(yyy), decltype(yyy1)>::value, "");
    static_assert(yyy == yyy1);
#if 1
    constexpr auto d = m / mm; // XXX radian !??!?!?!?!?!?!?!?!?!?!?!
    //using DDD = decltype(d);
    //static_assert(std::is_same_v<DDD::dimensions_type, units::kinds::Dimensionless<units::kinds::Length>>, "");
    constexpr double d_count = d.count();
    //constexpr double d_value = d.value();
#endif

#if 0
    constexpr auto fq1 = Hertzs{1 / 1_s};
    using FQ1 = decltype(fq1);
    static_assert(std::is_same_v<FQ1::dimension_type, units::kinds::Frequency>, "");
#endif

    ////constexpr Rainfalls xxxxxxxxx(1_cbm / 2_sqm);
    //////constexpr Rainfalls yyyyyyyyy = xxxxxxxxx + 2_cbm / 1_sqm;

    //////constexpr auto m001 = 1_cbm / 1_sqm;
    //////constexpr auto m011 = 1_cbm / 1_sqm;
    ////////constexpr auto m111 = m001 + m011;
    ////////constexpr auto m002 = m001 / 1_h;
    //////constexpr auto m002 = Rainfalls{m001} / 1_h;
    ////////constexpr auto m003 = 1_m + m001;
    ////////constexpr auto m003 = m001 + 1_m;
    ////////constexpr auto m003 = Rainfalls{m001} + 1_m;
    //////constexpr auto m004 = Metres{m001} + 1_m;
    ////////Incomplet<decltype(m001)>{};
    ////////constexpr auto m002 = m001 + 1_m;

    //constexpr auto vel2 = 1_m * fq1;
    //using VEL2 = decltype(vel2);
    //using VEL2dim = VEL2::dimension_type;
    //static_assert(std::is_same_v<VEL2::dimension_type, units::Velocity>, "");
#if 0
    constexpr auto d2 = d + 2.0;
    constexpr double d2_count = d2.count();
    constexpr double d2_value = d2.value();
    constexpr auto d3 = d * 2.0;
    constexpr double d3_count = d3.count();
    constexpr double d3_value = d3.value();
    constexpr bool cmp1 = d2 < d3;
    static_assert(cmp1);
#endif
    static_assert(1_cm != 1_m);
    static_assert(100_cm == 1_m);
    static_assert(1_m == 100_cm);
    static_assert(1_dm < 1_m);
    static_assert(1_m < 101_cm);
    static_assert(1_m > 99_cm);
    static_assert(101_cm > 1_m);
    static_assert(99_cm < 1_m);
#if 1
    constexpr bool cmp2 = 1_m >= 1_yd;
    static_assert(cmp2);
    static_assert(1_m != 1_yd);
    static_assert(1_m > 1_yd);
    static_assert(1_m >= 1_yd);
    static_assert(1_yd < 1_m);
    static_assert(1_yd <= 1_m);
    //static_assert(Min(1_yd, 1_m) == 1_yd);
    ////constexpr Metres min1 = Min(1_yd, 1_m);
    ////constexpr Yards min2 = Min(1_yd, 1_m);
    ////constexpr Millimetres min3 = Min(1_yd, 1_m);
    //constexpr Micrometres min4 = Min(1_yd, 1_m);
#endif
    //constexpr auto zzz = mm.convert_to(Millimetres::unit);
    //constexpr auto zxz = mm.convert_to(units::Millimeter{});
    constexpr auto aaa = mm.convert_to(Millimetres{});
    //constexpr auto bbb = mm.convert_to(Millimetres::unit_type{});
    //constexpr Metres x00 = m.convert_to(Millimetres{}) + 1_mm;
    //constexpr double x01 = 1001.0 * (1.0 / 1000.0);
    //constexpr double x02 = 1.0 / 1000.0;
    //using CM1 = CommonMultiplier< Millimetres::multiplier_type, Metres::multiplier_type >;
    //constexpr auto CM1_num = CM1::num;
    //constexpr auto CM2_den = CM1::den;
    //constexpr auto ddd01 = m / mm;
    //constexpr auto ddd02 = ddd01 * 3.4;
    //constexpr auto ddd03 = ddd02 + 5.6;
    //constexpr auto ddd04 = ddd03 / 3.4;
    //constexpr auto ddd05 = ddd04 - 5.6;
    //constexpr auto ddd06 = 3.4 * ddd05;
    //constexpr auto ddd07 = 5.6 + ddd06;
    //constexpr auto ddd08 = 3.4 / ddd07;
    //constexpr auto ddd09 = 5.6 - ddd08;

    //auto abc = 1_m;
    //abc += 1_cm;
    //abc += 1_mm;
    auto abc = 1_mm;
    abc += 1_cm;
    abc += 1_m;
    abc *= 2;
    auto def = 1_m / 1_m;
    //def += 2.0;
    def *= 2.0;
    def += 2_m / 2_m;
    //constexpr auto ghi = 2_m / 1_m + 1_m / 2_m;
    //constexpr auto ghi = 2_m + 50_cm;
    constexpr auto ghi = 1_m / 1_m + 1_cm / 1_cm;

    constexpr auto rad0 = Radians{1_rad / 2_rad};
//#if UNITS_HAS_NO_KINDS()
//    constexpr Radians phi0a = 1_rad / 2_rad;
//#endif
    constexpr auto rad1 = 1_rad + 2_rad /*+ rad0*/;
    //static_assert( units::kinds::IsDimensionless<Radians::dimension_type>::value, "");
    constexpr auto rad2 = Radians{1_m / 1_m};
    constexpr auto rad3 = rad1 + rad2;
    //constexpr auto rad4 = rad1 + 1_m / 1_m;
//#if !UNITS_IGNORE_KIND()
//    static_assert(!Compiles<Add, Radians, decltype(1_m / 1_m)>, "");
//#endif
    //constexpr auto rad5 = rad3 + 1_sr;
    //constexpr auto rad6 = rad1 + 1_m;
    constexpr auto deg1 = Degrees{1_rad};
    constexpr auto rad9 = Radians{1_deg};
    constexpr auto rad8 = Radians{1_gon};

    constexpr auto omega1 = 1_rad / 1_s;
    constexpr auto omega2 = 1 / 1_s;
    //constexpr auto omega3 = omega1 + omega2;

    constexpr auto lt1 = 1_m;
    constexpr auto lt2 = 1_s;
    constexpr auto lt3 = lt1 * lt2;
    constexpr auto lt4 = lt3 + 1_km * 1_h;

    ////constexpr auto rainfall = 1234_ml / 1_sqm; // rainfall ? should be m^3 / m^2 (not simplified) ?!?!?!
    ////constexpr auto rrr1 = rainfall / 1_h;
    //////constexpr auto rrr2 = 2.34_mps + rrr1; // Velocity + Rainfall / Hour ?!?!?!?!?!??!??!??!?!

    // constexpr auto E1 = 1_J;
    // constexpr auto E2 = 1_Nm;
    // constexpr auto E3 = Joules{E2};
    // //constexpr Joules E4 = 1_N * 1_m;
    // constexpr Joules E5 = Joules{1_N * 1_m};

#if 0
    constexpr double DDD = 1.234567890123456789e+13;
    const auto f0 = 1.0 / DDD;
    const auto f1 = 1.0 / f0;
    const auto f = DDD / Seconds{1};
    const auto ffff = 1_s / DDD;
    const auto f2 = 1 / (1_s / DDD);

    ////constexpr auto v1111 = Miles{12} / Hours{1};
    ////constexpr auto v2222 = 12 / Hours{1} * Miles{1};
    ////constexpr auto v3333 = quantity_cast<MetresPerSecond>(v2222);
    ////static_assert(v1111 == v2222);
    ////static_assert(v1111 == v3333);
    ////static_assert(v3333 == v1111);
    ////static_assert(std::is_same<decltype(v1111), decltype(v2222)>::value, "");
    //constexpr Quantity v1 = v1111;
    //constexpr Quantity v2 = quantity_cast<MetresPerSecond>(v1); // = v1.convert_to(metres / seconds);

    //constexpr MmPerSec yyy = v3;

    //constexpr auto zzzzz = Metres{1} + Seconds{2};

    ////constexpr Metres mm = 1_m;
    //////constexpr double yyy { 9832_mm };
    //////static_assert(yyy == 9832.0);
    //////constexpr double xxx { 1_m / 1_mm };
    //////static_assert(xxx == 1000.0);
    //////constexpr auto m2 = mm + (2 /** (kilometres * milliseconds / kilometres / milliseconds)*/) * Yards{1};
    ////constexpr auto m2 = 1_m + 2 * Yards{1}; // (2 * (1_km * 1_ms / 1_km / 1_ms)) * Yards{1};
    ////constexpr auto m3 = Metres{1} - mm;
    ////constexpr Yards yr = Metres{1};
    ////constexpr Metres m = Yards{1};
    //////constexpr auto msq = m * m + Yards{1};
#endif

    constexpr auto Q1 = 1_m / 1_m;
    constexpr auto Q2 = 1_s / 1_s;
//#if UNITS_IGNORE_KIND()
//    constexpr auto Q3a = Q1 + Q2;
//    constexpr auto Q3b = (1_m * 1_s + 1_s * 1_m) / (1_m * 1_s);
//    static_assert(std::is_same<decltype(Q3a), decltype(Q3b)>::value, "");
//    constexpr auto Q3c = (1_m * 1_s + 1_s * 1_m) / 1_m / 1_s;
//    static_assert(std::is_same<decltype(Q3a), decltype(Q3c)>::value, "");
//#endif
    constexpr auto Q3 = (1_m * 1_s + 1_m * 1_s) / (1_m * 1_s);
    //Incomplet<decltype(Q3)>{};
    //using XXXXX = decltype(Q3);

    constexpr auto area1 = 1_m * 1_m;
    constexpr auto area2 = 1_cm * 1_cm;
    constexpr auto area3 = 1_mi * 1_yd;
    constexpr auto area4 = area1 + area2 + area3;

    constexpr auto ouch1 = 1_rad + Radians{1_rad / 1_rad};

    //constexpr auto ouch2 = 1_m / 1_s + 1_s / 1_m;

    constexpr auto ouch2 = (1_m / 1_s) / (1_km / 1_h);
    //Incomplet<decltype(ouch2)>{};

    //constexpr auto ouch2 = 1_rad + 2.0;
    //constexpr auto ouch2a = 1_rad / 1_rad + 2.0;
    //constexpr auto ouch3 = 1_m + 2.0;
    //constexpr auto ouch4 = 1+ ((1_m/1_s) * (1_m/1_s)) / ((1_km/1_h) * (1_km/1_h));

    //constexpr auto speed1 = avg_speed(220_km, 2_h);
    //constexpr auto speed2 = avg_speed(140_mi, 2_h);

    //////constexpr auto ouch17 = KilometresPerHour{1_c} + MetresPerSecond{2_m / 2_s};
    //constexpr auto ouch18 = MetresPerSecond{0.25_c + 0.5_c} + 1_km / 1_h;
    //constexpr MetresPerSecond mps = 1_km / 1_h;
    //////constexpr MetresPerSecond vel01 = 1_c;
    //////constexpr SpeedOfLight vel02 = SpeedOfLight{299'792'458_m / 1_s};
    //////constexpr auto vel03 = 1_m / 1_s;

    //constexpr Micrometres imp01 = 1_in;
    //constexpr Millimetres imp02 = 1_in;
    //////constexpr auto vel01 = 1_m / 1_s + 1_km / 1_h;
    //////constexpr auto vel02 = kind_cast<kinds::Velocity>(1_m / 1_s) + 1_miph;

    auto xxxxx = 1_m / 1_s;
    auto yyyyy = MetresPerSecond{2_km / 1_h};
    xxxxx += kind_cast<DivKinds<kinds::Length, kinds::Time>>(yyyyy);
    xxxxx += 1_km / 1_s;
    //xxxxx += 1_km / 1_h;

    return 0;
}
#endif
