#include "../src/Units.h"

#include <stdio.h>
#include <chrono>

template <int L = 0, int M = 0>
struct Exp {};

template <typename D, typename E>
struct Dim {
    using dimension = D;
    using exponents = E;
};

struct Length : Dim<Length, Exp<1, 0>> {};
struct Mass : Dim<Mass, Exp<0, 1>> {};

template <typename...> struct Incomplet;

namespace units {
    namespace dimensions {
        struct Rainfall : Dimension<Rainfall, decltype(Volume::exponents{} / Area::exponents{})> {};

        //constexpr auto operator/(Volume, Area) { return Rainfall{}; }
    }

    using Rainfall  = Unit< Ratio<1, 1>, dimensions::Rainfall >;
    using Rainfalls = Quantity< Rainfall >;
}

namespace units {
    namespace dimensions {
        struct LengthTime : Dimension<LengthTime, decltype(Length::exponents{} * Time::exponents{})> {};

        constexpr auto operator*(Time, Length) { return LengthTime{}; }
        constexpr auto operator*(Length, Time) { return LengthTime{}; }
        constexpr auto operator/(LengthTime, Time) { return Length{}; }
        constexpr auto operator/(LengthTime, Length) { return Time{}; }
    }
}

namespace units {
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

int main()
{
    using namespace units;
    using namespace units::literals;

    constexpr auto hhhh = 1_h + 1_s;
    constexpr Meters m = 1_m;
    constexpr Millimeters mm = 1_mm;
    constexpr auto mmsq = mm * mm;
    constexpr auto mm01 = mm.count();
    constexpr auto mm02 = mm.value();
    constexpr auto m2 = 1_m + 1_cm + 1_mm + 1_nm;
    constexpr auto one_m = 1_m;
    constexpr auto two_yd = 2_yd;
    constexpr auto xxx = 1_m + 2_yd + 1_mm;
    constexpr auto xxx2 = 2_yd + 1_m + 1_mm;
    using TTT = decltype(xxx2);
    constexpr auto ttt = TTT::conversion;
    constexpr auto tttNum = TTT::conversion_type::num;
    constexpr auto tttDen = TTT::conversion_type::den;
    constexpr double xxx_value = xxx.value();
    constexpr auto yyy = quantity_cast<Millimeters>(xxx);
    constexpr auto yyy1 = Millimeters{xxx};
    constexpr Millimeters yyy2 = xxx;
    constexpr Millimeters yyy3 = 1_m;
    static_assert(std::is_same_v<decltype(yyy), decltype(yyy1)>);
    static_assert(yyy == yyy1);
#if 1
    constexpr auto d = m / mm; // XXX radian !??!?!?!?!?!?!?!?!?!?!?!
    //using DDD = decltype(d);
    //static_assert(std::is_same_v<DDD::dimensions_type, units::dimensions::Dimensionless<units::dimensions::Length>>, "");
    constexpr double d_count = d.count();
    constexpr double d_value = d.value();
#endif

    constexpr auto fq1 = 1 / 1_s;
    using FQ1 = decltype(fq1);
    static_assert(std::is_same_v<FQ1::dimension_type, units::dimensions::Frequency>, "");

    constexpr Rainfalls xxxxxxxxx(1_cbm / 2_sqm);
    //constexpr Rainfalls yyyyyyyyy = xxxxxxxxx + 2_cbm / 1_sqm;

    constexpr auto m001 = 1_cbm / 1_sqm;
    constexpr auto m011 = 1_cbm / 1_sqm;
    //constexpr auto m111 = m001 + m011;
    //constexpr auto m002 = m001 / 1_h;
    constexpr auto m002 = Rainfalls{m001} / 1_h;
    //constexpr auto m003 = 1_m + m001;
    //constexpr auto m003 = m001 + 1_m;
    //constexpr auto m003 = Rainfalls{m001} + 1_m;
    constexpr auto m004 = Meters{m001} + 1_m;
    //Incomplet<decltype(m001)>{};
    //constexpr auto m002 = m001 + 1_m;

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
    constexpr bool cmp2 = 1_m >= 1_yd;
    static_assert(cmp2);
    static_assert(1_m == 1000_mm);
    static_assert(1_m != 1_yd);
    static_assert(1_m > 1_yd);
    static_assert(1_m >= 1_yd);
    static_assert(1_yd < 1_m);
    static_assert(1_yd <= 1_m);
    constexpr auto zzz = mm.convert_to(Millimeters::unit);
    constexpr auto zxz = mm.convert_to(Millimeter{});
    constexpr auto aaa = mm.convert_to(Millimeters{});
    constexpr auto bbb = mm.convert_to(Millimeters::unit_type{});
    constexpr Meters x00 = m.convert_to(Millimeters{}) + 1_mm;
    constexpr double x01 = 1001.0 * (1.0 / 1000.0);
    constexpr double x02 = 1.0 / 1000.0;
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

    auto abc = 1_m;
    abc += 1_cm;
    abc += 1_mm;
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
    //static_assert( units::dimensions::IsDimensionless<Radians::dimension_type>::value, "");
    constexpr auto rad2 = Radians{1_m / 1_m};
    constexpr auto rad3 = rad1 + rad2;
    //constexpr auto rad4 = rad1 + 1_m / 1_m;
    //constexpr auto rad5 = rad3 + 1_sr;
    //constexpr auto rad6 = rad1 + 1_m;

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
    constexpr auto Q3 = (1_m * 1_s + 1_s * 1_m) / (1_m * 1_s);
    //Incomplet<decltype(Q3)>{};
    using XXXXX = decltype(Q3);

    return 0;
}
