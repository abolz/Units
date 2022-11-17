#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

static void test()
{
    using Width  = TaggedQuantity<Millimeters, struct _width>;
    using Height = TaggedQuantity<Centimeters, struct _height>;
    using Area   = decltype(Width{} * Height{});

    Width  w(1.0);
    Width  w2(1.0_cm);
    Width  w3(w2);
    Height h = convert_to<Height>(2.0_cm);

#if 0
    w = h; // will not compile
    // const auto ignored = w + h; // will not compile
#endif

    Area area = w * h; // (some kind of "area")
    //area = w * h * w;
    w = area / h; // works
    h = area / w; // works
    const auto value = count_as<SquareCentimeters>(area);

    const auto a = convert_to<SquareMillimeters>(area); // explicit cast works
}

static void test2()
{
    // static_assert( !std::is_constructible_v<SquareCentimetersPerMeter, Kilonewtons> );

    constexpr SquareCentimetersPerMeter as(SquareCentimeters(2.5) / Meters(1));
    //constexpr SquareCentimetersPerMeter as2 = SquareCentimeters(2.5) / Meters(1);
    //constexpr SquareCentimetersPerMeter as2 = SquareMeters(2.5) / Meters(1);
    //constexpr SquareCentimetersPerMeter as2 = SquareMeters(2.5) / Meters(1);
    constexpr SquareCentimetersPerMeter as2 = SquareCentimeters(2.5) / Meters(1);
    constexpr SquareCentimetersPerMeter as3(2.5);
    static_assert( as2 == as3 );
    //constexpr SquareCentimetersPerMeter as(Kilonewtons(2.5));
    // const auto xxx = as * Meters(1);
    // as = SquareCentimetersPerMeter( SquareCentimeters(5) / Meters(2) );

    //constexpr auto xxx = as + Kilonewtons(1);

    constexpr auto XXX = as2 + SquareCentimetersPerMeter(1);
    constexpr auto YYY = as2 + SquareCentimeters(1) / Meters(1);
    //constexpr auto ZZZ = as2 + 1.0;

    // constexpr Millimeters yyy = as;
    // constexpr Millimeters yyy(as);
    constexpr Millimeters yyy(as.untagged());
}

//static void test2()
//{
//    Watts w;
//    Vars v;
//    //Watts w0 = v;
//    //Vars v0 = w;
//    Watts w1(v);
//    Vars v1(w);
//}

//int main()
//{
//    test();
//}
