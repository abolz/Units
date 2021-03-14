#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

#include <limits>

using namespace uom;
using namespace uom::literals;

// f(x) = sin(x)
static auto f(Radians x)
{
    return sin(x);
}

// f'(x) = cos(x)
static auto f_prime(Radians x)
{
    // constexpr decltype(x) h(1.0e-10);
    // return (f(x + h) - f(x - h)) / (2 * h);

    return cos(x) / 1_rad;
}

// x_{n+1} = x_n - f(x_n) / f'(x_n)
static auto newton_step(Radians x)
{
    return x - f(x) / f_prime(x);
}

#if 0
int main()
{
    const auto s1 = newton_step(Radians(0.4_rev));
    const auto s2 = newton_step(s1);
    const auto s3 = newton_step(s2);
    const auto s4 = newton_step(s3);
    const auto s5 = newton_step(s4);
    const auto s6 = newton_step(s5);
    const auto s7 = newton_step(s6);
    const auto s8 = newton_step(s7);
    const auto s9 = newton_step(s8);

    return 0;
}
#endif
