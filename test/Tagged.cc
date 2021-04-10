#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

static void test()
{
    using Width  = TaggedQuantity<Millimeters, struct _width>;
    using Height = TaggedQuantity<Centimeters, struct _height>;

    Width  w(1.0);
    Height h = convert_to<Height>(2.0_cm);

#if 0
    w = h; // will not compile
    w + h; // will not compile
#endif

    const auto area = w * h; // (some kind of "area")
    w = area / h; // works
    h = area / w; // works
    const auto value = count_as<SquareCentimeters>(area);

    const auto a = convert_to<SquareMillimeters>(area); // explicit cast works
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
