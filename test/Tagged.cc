#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

static void test()
{
    using Width  = Tagged<Millimetres, struct _width>;
    using Height = Tagged<Centimetres, struct _height>;

    Width  w(1.0);
    Height h(2.0_cm);

#if 0
    w = h; // will not compile
    w + h; // will not compile
#endif

    const auto area = w * h; // (some kind of "area")
    w = area / h; // works
    h = area / w; // works
    const auto value = count<SquareCentimetres>(area);

    const auto a = SquareMillimetres(area); // explicit cast works
}

//int main()
//{
//    test();
//}
