#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

#include <limits>

using namespace uom;
using namespace uom::literals;

// f(x) = x^2
static constexpr auto f(Meters x)
{
    return x * x - 2.0_m2;
}

// f'(x) = 2x
static constexpr auto f_prime(Meters x)
{
    // constexpr decltype(x) h(1.0e-10);
    // return (f(x + h) - f(x - h)) / (2 * h);

    return 2 * x;
}

// x_{n+1} = x_n - f(x_n) / f'(x_n)
static constexpr auto newton_step(Meters x)
{
    return x - f(x) / f_prime(x);
}

static constexpr void testNewton()
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

static constexpr auto integrate_step(Meters x, Meters h, int n) -> Meters
{
    if (n <= 0)
        return 0_m;
    else
        return f_prime(x) + integrate_step(x + h, h, n - 1);
}

static constexpr auto integrate(Meters a, Meters b, int n)
{
    const auto h = (b - a) / n;

    return f_prime(a) * (h / 2)
         + integrate_step(a + h, h, n - 1) * h
         + f_prime(b) * (h / 2);
}

static constexpr void testIntegrate()
{
    //constexpr SquareMeters i1 = integrate(0_m, 2_m, 1);
    constexpr SquareMeters i2 = integrate(0_m, 2_m, 2);
    static_assert(i2._count_internal() == 4.0);

    //constexpr SquareMeters i3 = integrate(1_m, 2_m, 3);
    constexpr SquareMeters i4 = integrate(1_m, 2_m, 4);
    static_assert(i4._count_internal() == 3.0);

    //constexpr SquareMeters i5 = integrate(2_m, 3_m, 5);
    //constexpr SquareMeters i6 = integrate(2_m, 3_m, 6);
    //constexpr SquareMeters i7 = integrate(2_m, 3_m, 7);
    //constexpr SquareMeters i8 = integrate(2_m, 3_m, 8);
    constexpr SquareMeters i9 = integrate(2_m, 3_m, 9);
    static_assert(i9._count_internal() == 5.0);
}
