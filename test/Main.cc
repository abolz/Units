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
    constexpr Meters<> m1 = Millimeters<>{1};
    constexpr Millimeters<> m2 = Meters<>{1};
    constexpr auto m3 = m1 + m2;
    constexpr auto m4 = 1_m + 1_cm + 1_mm;
    //constexpr auto m5 = 1 / (1 / 2.98873687628376_m);
    //static_assert(m5 == 2.98873687628376_m, "");

    constexpr auto m6 = 1_mi + 1_yd + 1_ft + 1_in;
    constexpr auto m7 = Micrometers{m6};
    constexpr auto m8 = 1_mi + 1_yd;
    constexpr auto m9 = 1_yd + 1_ft;
    constexpr auto m10 = 1_ft + 1_in;
    constexpr auto m11 = 1_in + 1_cm;
    constexpr auto cm11 = Centimeters{m11};
    constexpr auto cm11x = 1.0 + 2.54;
    constexpr auto m12 = Centimeters{1_in};
    constexpr Centimeters m13 = 1_in;
    constexpr auto m14 = 1_in;
    //using InchType = units::Inch::conversion;
    //Incomplet<InchType>{};
    //Incomplet<decltype(m11)::conversion>{};
    constexpr auto m15 = Centimeters{Meters{Yards{Inches{1_cm}}}};
    //Incomplet<decltype(m15)::conversion>{};
    constexpr auto m16 = Centimeters{Inches{1_cm}};
    //constexpr auto m17 = Centimeters{Yards{Inches{1_cm}}};
    constexpr auto m18 = 3.2347687234_mi * 1_cm * 2.2347326487236487234_in / 1e+2_yd * 1e+2_yd / 2.2347326487236487234_in / 3.2347687234_mi;
    constexpr auto m18cm = Centimeters{m18};
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

#if 0
namespace sc {
    static_assert(Root(1, 1) == 1, "");
    static_assert(Root(1, 2) == 1, "");
    static_assert(Root(1, 10) == 1, "");
    static_assert(Root(1, INT64_MAX) == 1, "");
    static_assert(Root(1, INT64_MAX / 2) == 1, "");
    static_assert(Root(2, INT64_MAX) == 1, "");
    static_assert(Root(2, INT64_MAX / 2) == 1, "");
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
    static_assert(Root(INT64_MAX, 17) == 13, "");
    static_assert(Root(INT64_MAX, INT64_MAX) == 1, "");
    static_assert(Root(INT64_MAX, 2) == 3037000499, "");
    static_assert(Root(INT64_MAX, 3) == 2097151, "");
    static_assert(Root(INT64_MAX, 4) == 55108, "");
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

static void test001()
{
    constexpr auto x1 = 1_mps;
    constexpr auto y1 = kind_cast<kinds::Velocity>(1_m / 1_s);
    static_assert(x1 == y1);

    constexpr auto x2 = 1_mps;
    constexpr auto y2 = 1_m / 1_s;
    //static_assert(x2 == y2);
}

#if UNITS_HAS_ANY()
static constexpr int takesLength(Centimeters)
{
    return 0;
}

static void testFlatten()
{
    constexpr auto len1 = 1_m;
    constexpr auto dur1 = 1_s;
    //constexpr auto xxx1 = len1 * dur1 / dur1;
    constexpr auto xxx1 = len1 / dur1 * dur1;
    //constexpr auto res1 = takesLength(xxx1);
    constexpr auto res2 = takesLength(flatten(xxx1));
    //constexpr auto res3 = takesLength(as_any(xxx1));
    //constexpr auto res4 = takesLength(as_length(xxx1));
    //constexpr auto res5 = takesLength(as_length(1_s));
    //constexpr auto xxx2 = kind_cast<kinds::Length>(xxx1);
    //constexpr auto res6 = takesLength(xxx2);
    //constexpr auto res6 = takesLength(as_any_quantity(1_m)); // any_kind_cast(1_s));
    constexpr auto res7 = takesLength(1_m);
}
#endif

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
    constexpr auto sum0 = 1_mm * 1_h + 2_km * 1_h;
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

#if 1
static void test1()
{
    //constexpr auto x = 1_m / 1_m;
    //constexpr auto y = 1_s / 1_s;
    //constexpr auto z = x + y;

    //constexpr Meters zzz = 1_s;
    constexpr Centimeters zzz = 1_m;
    constexpr Centimeters zzz2{1_m};
    //constexpr Centimeters zzz3{1_s};
    //constexpr Meters z2 = Seconds{1}.convert_to(Meters{});

    auto vel1 = 1_mps;
    auto vel1a = 1_m * (1 / 1_s);
    auto vel1b = vel1 + kind_cast<kinds::Velocity>(vel1a);
    auto vel2 = 1_kmph;
    //vel1 = vel2;
    vel1 = MetersPerSecond{vel2};

    //auto len1 = 1_m;
    //len1 = 1_mm;
    auto len2 = 1_mm;
    len2 = 1_m;

    //constexpr Quantity m = 1_m;
    //constexpr Quantity s = 1_s;
    //constexpr Quantity z = m + s;

    auto vel01 = 1_kmph;
    //KilometersPerHour vel02 = 1_mps;
    KilometersPerHour vel03 = KilometersPerHour{1_mps};

    constexpr auto phi0 = 1_gon + 1_deg;
    //Incomplet<decltype(phi0)>{};
    constexpr auto phi1 = Degrees{phi0};
    //constexpr auto phi2 = 1 / 1_gon + 1 / 1_deg;
    constexpr auto phi3 = 1 / 1_gon;
    constexpr auto val3 = phi3.value();
    constexpr auto phi4 = 1 / 1_deg;
    constexpr auto val4 = phi4.value();
    constexpr auto phi5 = phi3 + phi4;
    constexpr auto val5 = phi5.value();
    //Incomplet<decltype(phi5)>{};

    //constexpr auto phi6 = 1_deg + 1_rad;
    //Incomplet<decltype(phi6)>{};

    //constexpr auto phi7 = 1_deg + 1_rev;
    //Incomplet<decltype(phi7)>{};
    //constexpr auto phi8 = 1_deg * 1_rev;
    //Incomplet<decltype(phi8)>{};
    //constexpr auto phi9 = phi7 / phi8;
    //Incomplet<decltype(phi9)>{};
    //constexpr auto phi10 = 1 / 1_deg + 1 / 1_rev;
    //Incomplet<decltype(phi10)>{};

    constexpr auto phi7 = 1_m + 1_km;
    //Incomplet<decltype(phi7)>{};
    constexpr auto phi8 = 1_m * 1_km;
    //Incomplet<decltype(phi8)>{};
    constexpr auto phi9 = phi7 / phi8;
    //Incomplet<decltype(phi9)>{};
    constexpr auto phi10 = 1 / 1_m + 1 / 1_km;
    //Incomplet<decltype(phi10)>{};
    //constexpr auto phi11 = phi9 - phi10;
    //static_assert(phi9 == phi10, "");
    using Phi9 = std::remove_const_t<decltype(phi9)>;
    constexpr auto phi12 = phi9 - Phi9(phi10);
    //Incomplet<decltype(phi12)>{};
#if UNITS_HAS_ANY()
    constexpr auto phi13 = flatten(phi9) - flatten(phi10);
    //Incomplet<decltype(phi13)>{};

    constexpr Phi9 phi14 = flatten(phi10);

    //constexpr auto phi15 = flatten(phi9) - phi10;
    //Incomplet<decltype(phi15)>{};

    constexpr auto phi16 = flatten(flatten(phi9) * flatten(phi10));
    //Incomplet<decltype(phi16)>{};
#endif
}
#endif

int main()
{
    constexpr auto hhhh = 1_h + 1_s;
    constexpr Meters m = 1_m;
    constexpr Meters m2222 = Meters{1_mm};
    constexpr Millimeters mm = 1_mm;
    constexpr Millimeters mmmm22222 = 1_m;
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
    constexpr auto yyy = quantity_cast<Millimeters>(xxx);
    constexpr auto yyy1 = Millimeters{xxx};
    //constexpr Millimeters yyy2 = xxx;
    constexpr Millimeters yyy3 = 1_m;
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
    //////constexpr auto m004 = Meters{m001} + 1_m;
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
    ////constexpr Meters min1 = Min(1_yd, 1_m);
    ////constexpr Yards min2 = Min(1_yd, 1_m);
    ////constexpr Millimeters min3 = Min(1_yd, 1_m);
    //constexpr Micrometers min4 = Min(1_yd, 1_m);
#endif
    //constexpr auto zzz = mm.convert_to(Millimeters::unit);
    //constexpr auto zxz = mm.convert_to(units::Millimeter{});
    constexpr auto aaa = mm.convert_to(Millimeters{});
    //constexpr auto bbb = mm.convert_to(Millimeters::unit_type{});
    //constexpr Meters x00 = m.convert_to(Millimeters{}) + 1_mm;
    //constexpr double x01 = 1001.0 * (1.0 / 1000.0);
    //constexpr double x02 = 1.0 / 1000.0;
    //using CM1 = CommonMultiplier< Millimeters::multiplier_type, Meters::multiplier_type >;
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
    constexpr auto rad1 = 1_rad + 2_rad /*+ rad0*/;
    //static_assert( units::kinds::IsDimensionless<Radians::dimension_type>::value, "");
    constexpr auto rad2 = Radians{1_m / 1_m};
    constexpr auto rad3 = rad1 + rad2;
    //constexpr auto rad4 = rad1 + 1_m / 1_m;
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
    ////constexpr auto v3333 = quantity_cast<MetersPerSecond>(v2222);
    ////static_assert(v1111 == v2222);
    ////static_assert(v1111 == v3333);
    ////static_assert(v3333 == v1111);
    ////static_assert(std::is_same<decltype(v1111), decltype(v2222)>::value, "");
    //constexpr Quantity v1 = v1111;
    //constexpr Quantity v2 = quantity_cast<MetersPerSecond>(v1); // = v1.convert_to(meters / seconds);

    //constexpr MmPerSec yyy = v3;

    //constexpr auto zzzzz = Meters{1} + Seconds{2};

    ////constexpr Meters mm = 1_m;
    //////constexpr double yyy { 9832_mm };
    //////static_assert(yyy == 9832.0);
    //////constexpr double xxx { 1_m / 1_mm };
    //////static_assert(xxx == 1000.0);
    //////constexpr auto m2 = mm + (2 /** (kilometers * milliseconds / kilometers / milliseconds)*/) * Yards{1};
    ////constexpr auto m2 = 1_m + 2 * Yards{1}; // (2 * (1_km * 1_ms / 1_km / 1_ms)) * Yards{1};
    ////constexpr auto m3 = Meters{1} - mm;
    ////constexpr Yards yr = Meters{1};
    ////constexpr Meters m = Yards{1};
    //////constexpr auto msq = m * m + Yards{1};
#endif

    constexpr auto Q1 = 1_m / 1_m;
    constexpr auto Q2 = 1_s / 1_s;
    //constexpr auto Q3 = Q1 + Q2;
    // constexpr auto Q3 = (1_m * 1_s + 1_s * 1_m) / (1_m * 1_s);
    constexpr auto Q3 = (1_m * 1_s + 1_m * 1_s) / (1_m * 1_s);
    //Incomplet<decltype(Q3)>{};
    //using XXXXX = decltype(Q3);

    constexpr auto area1 = 1_m * 1_m;
    constexpr auto area2 = 1_cm * 1_cm;
    constexpr auto area3 = area1 + area2 + 1_mi * 1_yd;

    constexpr auto ouch1 = 1_rad + Radians{1_rad / 1_rad};

    //constexpr auto ouch2 = 1_m / 1_s + 1_s / 1_m;

    constexpr auto ouch2 = (1_m / 1_s) / (1_km / 1_h);
    //Incomplet<decltype(ouch2)>{};

    //constexpr auto ouch2 = 1_rad + 2.0;
    //constexpr auto ouch3 = 1_m + 2.0;
    //constexpr auto ouch4 = 1+ ((1_m/1_s) * (1_m/1_s)) / ((1_km/1_h) * (1_km/1_h));

    //constexpr auto speed1 = avg_speed(220_km, 2_h);
    //constexpr auto speed2 = avg_speed(140_mi, 2_h);

    //////constexpr auto ouch17 = KilometersPerHour{1_c} + MetersPerSecond{2_m / 2_s};
    //constexpr auto ouch18 = MetersPerSecond{0.25_c + 0.5_c} + 1_km / 1_h;
    //constexpr MetersPerSecond mps = 1_km / 1_h;
    //////constexpr MetersPerSecond vel01 = 1_c;
    //////constexpr SpeedOfLight vel02 = SpeedOfLight{299'792'458_m / 1_s};
    //////constexpr auto vel03 = 1_m / 1_s;

    //constexpr Micrometers imp01 = 1_in;
    //constexpr Millimeters imp02 = 1_in;
    constexpr auto vel01 = 1_m / 1_s + 1_km / 1_h;
    constexpr auto vel02 = kind_cast<kinds::Velocity>(1_m / 1_s) + 1_miph;

    auto xxxxx = 1_m / 1_s;
    auto yyyyy = MetersPerSecond{2_km / 1_h};
    xxxxx += kind_cast<DivKinds<kinds::Length, kinds::Time>>(yyyyy);
    xxxxx += 1_km / 1_s;
    //xxxxx += 1_km / 1_h;

    return 0;
}
#endif
