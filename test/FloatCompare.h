#pragma once

namespace uom {

inline constexpr double kDefaultEps = 0x1p-48; // ~= 3.5 10^-15

constexpr double constexpr_min(double x, double y)
{
    return y < x ? y : x;
}

constexpr double constexpr_max(double x, double y)
{
    return y < x ? x : y;
}

constexpr double constexpr_abs(double x)
{
    return x < 0 ? -x : x;
}

constexpr bool almost_zero(double x, double absolute_eps)
{
    return constexpr_abs(x) <= absolute_eps;
}

// Answers the question:
//  "Does the smaller value lie within some neighborhood of the larger value?"
constexpr bool approximately_equal(double x, double y, double relative_eps = kDefaultEps)
{
    const double delta = constexpr_abs(x - y);
    const double abs_x = constexpr_abs(x);
    const double abs_y = constexpr_abs(y);

    return delta <= relative_eps * constexpr_max(abs_x, abs_y);
}

// Answers the question:
//  "Does the larger value lie within some neighborhood of the smaller value?"
constexpr bool essentially_equal(double x, double y, double relative_eps = kDefaultEps)
{
    const double delta = constexpr_abs(x - y);
    const double abs_x = constexpr_abs(x);
    const double abs_y = constexpr_abs(y);

    return delta <= relative_eps * constexpr_min(abs_x, abs_y);
}

constexpr bool almost_equal(double x, double y, double eps = kDefaultEps)
{
    // The absolute test fails when x and y are large,
    // the relative test fails when x and y are small.
    // The following test is equivalent to:
    //  |x-y| <= eps  OR  |x-y| <= eps * max(|x|, |y|)

    const double delta = constexpr_abs(x - y);
    const double abs_x = constexpr_abs(x);
    const double abs_y = constexpr_abs(y);

//  return delta <= constexpr_max(eps, eps * constexpr_max(abs_x, abs_y));
    return delta <= eps * constexpr_max(1.0, constexpr_max(abs_x, abs_y));
}

} // namespace uom
