// Copyright 2021 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "Unit.h"

#include <cmath>

namespace uom {

//==================================================================================================
//
//==================================================================================================

namespace impl
{
    constexpr bool MulOverflows(const int64_t x, const int64_t y) noexcept
    {
        UNITS_ASSERT(x >= 0);
        UNITS_ASSERT(y >= 0);

        return x > INT64_MAX / (y == 0 ? 1 : y);
    }

    // Computes y^n.
    // Does not check for overflow.
    constexpr int64_t Power(const int64_t y, const int64_t n) noexcept
    {
        UNITS_ASSERT(y >= 1);
        UNITS_ASSERT(n >= 0);

        if (n == 0)
            return 1;

        int64_t r = 1;
        int64_t a = y;
        int64_t k = n;

        while (k >= 2)
        {
            if (k % 2 != 0)
            {
                UNITS_ASSERT(!MulOverflows(r, a));
                r *= a;
            }
            UNITS_ASSERT(!MulOverflows(a, a));
            a *= a;
            k /= 2;
        }
        UNITS_ASSERT(!MulOverflows(r, a));
        r *= a;

        return r;
    }

    constexpr int Compare3Way(const int64_t x, const int64_t y) noexcept
    {
        if (x < y)
            return -1;
        if (x > y)
            return +1;

        return 0;
    }

    // Returns y^n <=> x (without overflow).
    constexpr int ComparePower(const int64_t y, const int64_t n, const int64_t x) noexcept
    {
        UNITS_ASSERT(y >= 1);
        UNITS_ASSERT(n >= 0);
        UNITS_ASSERT(x >= 0);

        if (n == 0) // r == 1
            return Compare3Way(1, x);

        int64_t r = 1;
        int64_t a = y;
        int64_t k = n;

        while (k >= 2)
        {
            if (k % 2 != 0)
            {
                if (MulOverflows(r, a))
                    return +1;
                r *= a;
            }
            if (MulOverflows(a, a))
                return +1;
            a *= a;
            k /= 2;
        }
        if (MulOverflows(r, a))
            return +1;
        r *= a;

        return Compare3Way(r, x);
    }

    struct RootResult
    {
        int64_t value;
        bool is_exact;
    };

    // Computes y = n-th root of x,
    // i.e. returns the largest y, such that y^n <= x
    constexpr RootResult Root(int64_t x, int64_t n) noexcept
    {
        UNITS_ASSERT(x >= 0);
        UNITS_ASSERT(n >= 1);

        if (x <= 1 || n <= 1)
            return {x, true};

        int64_t lo = 1;
        int64_t hi = 1 + x / n;
        // Bernoulli  ==>  x^(1/n) <= 1 + (x - 1)/n < 1 + x/n
        // Since n >= 2, hi will not overflow here.

        for (;;)
        {
            const auto y = lo + (hi - lo) / 2;
            const auto cmp = ComparePower(y, n, x);

            if (y == lo || cmp == 0)
                return {y, cmp == 0};

            if (cmp < 0)
                lo = y;
            else if (0 < cmp)
                hi = y;
        }
    }

    template <typename R, int64_t N>
    struct RationalPower
    {
        static constexpr auto pow_num = impl::Power(R::num, N);
        static constexpr auto pow_den = impl::Power(R::den, N);

        using ratio = Ratio<pow_num, pow_den>;
    };

    template <typename R, int64_t N>
    struct RationalRoot
    {
        static constexpr auto root_num = impl::Root(R::num, N);
        static constexpr auto root_den = impl::Root(R::den, N);

        using ratio = Ratio<root_num.value, root_den.value>;
        static constexpr bool is_exact = (root_num.is_exact && root_den.is_exact);
    };

} // namespace impl

//==================================================================================================
// abs
//==================================================================================================

template <typename U>
[[nodiscard]] inline Quantity<U> abs(Quantity<U> q) noexcept
{
    return Quantity<U>(std::abs(q.count_internal()));
}

template <typename U>
[[nodiscard]] inline Quantity<U> fabs(Quantity<U> q) noexcept
{
    return Quantity<U>(std::fabs(q.count_internal()));
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Q, typename Z>
[[nodiscard]] inline Absolute<Q, Z> abs(Absolute<Q, Z> q) noexcept
{
    return Absolute<Q, Z>(std::abs(q.count_internal()));
}

template <typename Q, typename Z>
[[nodiscard]] inline Absolute<Q, Z> fabs(Absolute<Q, Z> q) noexcept
{
    return Absolute<Q, Z>(std::fabs(q.count_internal()));
}

//==================================================================================================
// min/max
//==================================================================================================

template <typename U>
[[nodiscard]] inline constexpr Quantity<U> min(Quantity<U> x, Quantity<U> y) noexcept
{
    // std::min
    return y.count_internal() < x.count_internal() ? y : x;
}

template <typename U>
[[nodiscard]] inline constexpr Quantity<U> max(Quantity<U> x, Quantity<U> y) noexcept
{
    // std::max
    return y.count_internal() < x.count_internal() ? x : y;
}

template <typename U>
[[nodiscard]] inline Quantity<U> fmin(Quantity<U> x, Quantity<U> y) noexcept
{
    return Quantity<U>(std::fmin(x.count_internal(), y.count_internal()));
}

template <typename U>
[[nodiscard]] inline Quantity<U> fmax(Quantity<U> x, Quantity<U> y) noexcept
{
    return Quantity<U>(std::fmax(x.count_internal(), y.count_internal()));
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Q, typename Z>
[[nodiscard]] inline constexpr Absolute<Q, Z> min(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    // std::min
    return y.count_internal() < x.count_internal() ? y : x;
}

template <typename Q, typename Z>
[[nodiscard]] inline constexpr Absolute<Q, Z> max(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    // std::max
    return y.count_internal() < x.count_internal() ? x : y;
}

template <typename Q, typename Z>
[[nodiscard]] inline Absolute<Q, Z> fmin(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    return Absolute<Q, Z>(std::fmin(x.count_internal(), y.count_internal()));
}

template <typename Q, typename Z>
[[nodiscard]] inline Absolute<Q, Z> fmax(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    return Absolute<Q, Z>(std::fmax(x.count_internal(), y.count_internal()));
}

//==================================================================================================
// pow
//==================================================================================================

// q^N
template <int64_t N, typename U>
[[nodiscard]] inline auto pow(Quantity<U> q)
{
    using conversion    = typename U::conversion;   // Conversion<Ratio, exp>
    using kind          = typename U::kind;         // Kind<dimension, Tag>
    using dimension     = typename U::dimension;    // Ratio

    static_assert(N >= 1,
        "only positive exponents allowed");
    static_assert(std::is_same<typename kind::tag, kinds::Simple>::value,
        "operation not supported - instead of 'pow<N>(q)' you must use 'pow<N>(q.simplify())'");
    //
    // TODO:
    // Complex-tags?
    //

    using conversion_pow = impl::RationalPower<typename conversion::ratio, N>;
    using dimension_pow = impl::RationalPower<dimension, N>;

    using C = Conversion<typename conversion_pow::ratio, conversion::exp * N>;
    using D = typename dimension_pow::ratio;

    using K = Kind<D, kinds::Simple>;

    if constexpr (N == 1)
    {
        return Quantity<Unit<C, K>>(q.count_internal());
    }
    else if constexpr (N == 2)
    {
        return Quantity<Unit<C, K>>(q.count_internal() * q.count_internal());
    }
    else
    {
        return Quantity<Unit<C, K>>(std::pow(q.count_internal(), static_cast<double>(N)));
    }
}

// q^2
template <typename U>
[[nodiscard]] inline auto square(Quantity<U> q)
{
    return uom::pow<2>(q);
}

//==================================================================================================
// nth_root
//==================================================================================================

// q^(1/N)
template <int64_t N, typename U>
[[nodiscard]] inline auto nth_root(Quantity<U> q)
{
    using conversion    = typename U::conversion;   // Conversion<Ratio, exp>
    using kind          = typename U::kind;         // Kind<dimension, Tag>
    using dimension     = typename U::dimension;    // Ratio

    static_assert(N >= 1,
        "only positive exponents allowed");
    static_assert(std::is_same<typename kind::tag, kinds::Simple>::value,
        "operation not supported - instead of 'nth_root<N>(q)' you must use 'nth_root<N>(q.simplify())'");
    //
    // TODO:
    // Complex-tags?
    //

    using conversion_root = impl::RationalRoot<typename conversion::ratio, N>;
    using dimension_root = impl::RationalRoot<dimension, N>;

    static_assert(conversion_root::is_exact,
        "inexact roots are not supported");
    static_assert(conversion::exp % N == 0,
        "inexact roots are not supported");
    static_assert(dimension_root::is_exact,
        "rational roots are not supported");

    using C = Conversion<typename conversion_root::ratio, conversion::exp / N>;
    using D = typename dimension_root::ratio;

    using K = Kind<D, kinds::Simple>;

    if constexpr (N == 1)
    {
        return Quantity<Unit<C, K>>(q.count_internal());
    }
    else if constexpr (N == 2)
    {
        return Quantity<Unit<C, K>>(std::sqrt(q.count_internal()));
    }
    else if constexpr (N == 3)
    {
        return Quantity<Unit<C, K>>(std::cbrt(q.count_internal()));
    }
    else
    {
        return Quantity<Unit<C, K>>(std::pow(q.count_internal(), 1.0 / static_cast<double>(N)));
    }
}

// q^(1/2)
template <typename U>
[[nodiscard]] inline auto sqrt(Quantity<U> q)
{
    return uom::nth_root<2>(q);
}

// q^(1/3)
template <typename U>
[[nodiscard]] inline auto cbrt(Quantity<U> q)
{
    return uom::nth_root<3>(q);
}

//==================================================================================================
// hypot
//==================================================================================================

template <typename U>
[[nodiscard]] inline Quantity<U> hypot(Quantity<U> x, Quantity<U> y) noexcept
{
    return Quantity<U>(std::hypot(x.count_internal(), y.count_internal()));
}

template <typename U>
[[nodiscard]] inline Quantity<U> hypot(Quantity<U> x, Quantity<U> y, Quantity<U> z) noexcept
{
    return Quantity<U>(std::hypot(x.count_internal(), y.count_internal(), z.count_internal()));
}

//==================================================================================================
// exp
//==================================================================================================

[[nodiscard]] inline Dimensionless exp(Dimensionless q) noexcept
{
    return Dimensionless(std::exp(q.count_internal()));
}

[[nodiscard]] inline Dimensionless exp2(Dimensionless q) noexcept
{
    return Dimensionless(std::exp2(q.count_internal()));
}

[[nodiscard]] inline Dimensionless expm1(Dimensionless q) noexcept
{
    return Dimensionless(std::expm1(q.count_internal()));
}

[[nodiscard]] inline Dimensionless log(Dimensionless q) noexcept
{
    return Dimensionless(std::log(q.count_internal()));
}

[[nodiscard]] inline Dimensionless log10(Dimensionless q) noexcept
{
    return Dimensionless(std::log10(q.count_internal()));
}

[[nodiscard]] inline Dimensionless log1p(Dimensionless q) noexcept
{
    return Dimensionless(std::log1p(q.count_internal()));
}

[[nodiscard]] inline Dimensionless log2(Dimensionless q) noexcept
{
    return Dimensionless(std::log2(q.count_internal()));
}

//==================================================================================================
// Trigonometric
//==================================================================================================

#if 1
[[nodiscard]] inline Dimensionless sin(Radians x) noexcept
{
    return Dimensionless(std::sin(x.count_internal()));
}

[[nodiscard]] inline Dimensionless cos(Radians x) noexcept
{
    return Dimensionless(std::cos(x.count_internal()));
}

[[nodiscard]] inline Dimensionless tan(Radians x) noexcept
{
    return Dimensionless(std::tan(x.count_internal()));
}
#else
template <typename U>
[[nodiscard]] inline Dimensionless sin(Quantity<U> x) noexcept
{
    static_assert(typename U::dimension == dims::PlaneAngle,
        "operation not supported");

    return Dimensionless(std::sin(x.template count<Radians>()));
}

template <typename U>
[[nodiscard]] inline Dimensionless cos(Quantity<U> x) noexcept
{
    static_assert(typename U::dimension == dims::PlaneAngle,
        "operation not supported");

    return Dimensionless(std::cos(x.template count<Radians>()));
}

template <typename U>
[[nodiscard]] inline Dimensionless tan(Quantity<U> x) noexcept
{
    static_assert(typename U::dimension == dims::PlaneAngle,
        "operation not supported");

    return Dimensionless(std::tan(x.template count<Radians>()));
}
#endif

[[nodiscard]] inline Radians asin(Dimensionless x) noexcept
{
    return Radians(std::asin(x.count_internal()));
}

[[nodiscard]] inline Radians acos(Dimensionless x) noexcept
{
    return Radians(std::acos(x.count_internal()));
}

[[nodiscard]] inline Radians atan(Dimensionless x) noexcept
{
    return Radians(std::atan(x.count_internal()));
}

#if 1
template <typename C1, typename C2, typename K>
[[nodiscard]] inline Radians atan2(Quantity<Unit<C1, K>> y, Quantity<Unit<C2, K>> x) noexcept
{
    using Q = Quantity<Unit<CommonConversion<C1, C2>, kinds::Simple>>;
    return Radians(std::atan2(Q(y).count_internal(), Q(x).count_internal()));
}
#else
[[nodiscard]] inline Radians atan2(Dimensionless y, Dimensionless x) noexcept
{
    return Radians(std::atan2(y.count_internal(), x.count_internal()));
}
#endif

[[nodiscard]] inline Dimensionless sinh(Dimensionless x) noexcept
{
//  return (exp(x) - exp(-x)) / 2;
    return Dimensionless(std::sinh(x.count_internal()));
}

[[nodiscard]] inline Dimensionless cosh(Dimensionless x) noexcept
{
//  return (exp(x) + exp(-x)) / 2;
    return Dimensionless(std::cosh(x.count_internal()));
}

[[nodiscard]] inline Dimensionless tanh(Dimensionless x) noexcept
{
//  return sinh(x) / tanh(x);
    return Dimensionless(std::tanh(x.count_internal()));
}

[[nodiscard]] inline Dimensionless asinh(Dimensionless x) noexcept
{
    return Dimensionless(std::asinh(x.count_internal()));
}

[[nodiscard]] inline Dimensionless acosh(Dimensionless x) noexcept
{
    return Dimensionless(std::acosh(x.count_internal()));
}

[[nodiscard]] inline Dimensionless atanh(Dimensionless x) noexcept
{
    return Dimensionless(std::atanh(x.count_internal()));
}

//==================================================================================================
//
//==================================================================================================

namespace impl
{
    template <typename Scale, typename U, typename Fn>
    inline auto ToInt(Quantity<U> x, Fn fn) noexcept
    {
        if constexpr (IsRatio<Scale>)
        {
            using T = ScaledQuantity<Conversion<Scale>, Quantity<U>>;
            return convert_to<Quantity<U>>( T( fn( convert_to<T>(x).count_internal() ) ) );
        }
        else // if constexpr (IsQuantity<Scale>)
        {
            using T = Scale;
            return T( fn( convert_to<T>(x).count_internal() ) );
        }
    }

    template <typename Scale, typename Q, typename Z, typename Fn>
    inline auto ToInt(Absolute<Q, Z> x, Fn fn) noexcept
    {
        return Absolute<Q, Z>( ToInt<Scale>( Q(x.count_internal()) ).count_internal() );
    }
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Scale = Ratio<1>, typename U>
[[nodiscard]] inline auto ceil(Quantity<U> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const double t) { return std::ceil(t); } );
}

template <typename Scale = Ratio<1>, typename U>
[[nodiscard]] inline auto floor(Quantity<U> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const double t) { return std::floor(t); } );
}

template <typename Scale = Ratio<1>, typename U>
[[nodiscard]] inline auto trunc(Quantity<U> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const double t) { return std::trunc(t); } );
}

template <typename Scale = Ratio<1>, typename U>
[[nodiscard]] inline auto round(Quantity<U> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const double t) { return std::round(t); } );
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Scale = Ratio<1>, typename Q, typename Z>
[[nodiscard]] inline auto ceil(Absolute<Q, Z> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const double t) { return std::ceil(t); } );
}

template <typename Scale = Ratio<1>, typename Q, typename Z>
[[nodiscard]] inline auto floor(Absolute<Q, Z> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const double t) { return std::floor(t); } );
}

template <typename Scale = Ratio<1>, typename Q, typename Z>
[[nodiscard]] inline auto trunc(Absolute<Q, Z> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const double t) { return std::trunc(t); } );
}

template <typename Scale = Ratio<1>, typename Q, typename Z>
[[nodiscard]] inline auto round(Absolute<Q, Z> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const double t) { return std::round(t); } );
}

//==================================================================================================
//
//==================================================================================================

template <typename U, typename U2, std::enable_if_t<std::is_same_v<MulUnits<U, U>, U2>, int> = 0>
[[nodiscard]] inline auto fma(Quantity<U> x, Quantity<U> y, Quantity<U2> z)
{
    // The function signature is more restrictive than (x * y + z).
    // This is to avoid any implicit conversion taking place, which would defeat the purpose
    // the fma instruction.

    return Quantity<MulUnits<U, U>>(std::fma(x.count_internal(), y.count_internal(), z.count_internal()));
}

//==================================================================================================
// midpoint
//==================================================================================================

template <typename U>
[[nodiscard]] inline Quantity<U> midpoint(Quantity<U> x, Quantity<U> y) noexcept
{
    return Quantity<U>(x.count_internal() / 2 + y.count_internal() / 2);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Q, typename Z>
[[nodiscard]] inline Absolute<Q, Z> midpoint(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    return Absolute<Q, Z>(x.count_internal() / 2 + y.count_internal() / 2);
}

//==================================================================================================
// lerp
//==================================================================================================

template <typename U>
[[nodiscard]] inline Quantity<U> lerp(Quantity<U> x, Quantity<U> y, double t) noexcept
{
    return Quantity<U>((1 - t) * x.count_internal() + t * y.count_internal());
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Q, typename Z>
[[nodiscard]] inline Absolute<Q, Z> lerp(Absolute<Q, Z> x, Absolute<Q, Z> y, double t) noexcept
{
    return Absolute<Q, Z>((1 - t) * x.count_internal() + t * y.count_internal());
}

//==================================================================================================
//
//==================================================================================================

namespace impl
{
    inline constexpr double kPi    = 3.14159265358979323846264338328;
    inline constexpr double kTwoPi = 6.28318530717958647692528676656;

    inline double Clamp(double x, double lo, double hi)
    {
        if (x < lo)
            return lo;
        if (x > hi)
            return hi;

        return x;
    }

    inline double NormalizePositive(double x, double range) noexcept
    {
        // UNITS_ASSERT(std::isfinite(x));

        double w = x;
        if (w < 0 || w >= range)
        {
            w = std::fmod(w, range);
            // POST: -range < w < range

            if (w < 0)
                w += range;

            w = Clamp(w, 0.0, range);
            if (w == range)
                w = 0;

            UNITS_ASSERT(w >= 0);
            UNITS_ASSERT(w < range);
        }

        return w;
    }

    inline double NormalizeSymmetric(double x, double range) noexcept
    {
        // UNITS_ASSERT(std::isfinite(x));

        const double half = range / 2.0;

        double w = x;
        if (w < -half || w >= half)
        {
            w = std::fmod(x, range);
            // POST: -range < w < range

            if (w < -half)
                w += range;
            else if (w > half)
                w -= range;

            w = Clamp(w, -half, half);
            if (w == half)
                w = -half;

            UNITS_ASSERT(w >= -half);
            UNITS_ASSERT(w < half);
        }

        return w;
    }

} // namespace impl

// Reduce x to [0, 2pi)
inline Radians normalize_positive(Radians x) noexcept
{
    return Radians(impl::NormalizePositive(x.count_internal(), impl::kTwoPi));
}

// Reduce x to [-pi, pi)
inline Radians normalize_symmetric(Radians x) noexcept
{
    return Radians(impl::NormalizeSymmetric(x.count_internal(), impl::kTwoPi));
}

// Reduce x to [0, 360)
inline Degrees normalize_positive(Degrees x) noexcept
{
    return Degrees(impl::NormalizePositive(x.count_internal(), 360.0));
}

// Reduce x to [-180, 180)
inline Degrees normalize_symmetric(Degrees x) noexcept
{
    return Degrees(impl::NormalizeSymmetric(x.count_internal(), 360.0));
}

// Reduce x to [0, 400)
inline Gons normalize_positive(Gons x) noexcept
{
    return Gons(impl::NormalizePositive(x.count_internal(), 400.0));
}

// Reduce x to [-200, 200)
inline Gons normalize_symmetric(Gons x) noexcept
{
    return Gons(impl::NormalizeSymmetric(x.count_internal(), 400.0));
}

// Reduce x to [0, 1)
inline Revolutions normalize_positive(Revolutions x) noexcept
{
    return Revolutions(impl::NormalizePositive(x.count_internal(), 1.0));
}

// Reduce x to [-1/2, 1/2)
inline Revolutions normalize_symmetric(Revolutions x) noexcept
{
    return Revolutions(impl::NormalizeSymmetric(x.count_internal(), 1.0));
}

} // namespace uom
