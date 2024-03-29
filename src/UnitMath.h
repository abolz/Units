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
[[nodiscard]] constexpr Quantity<U> abs(Quantity<U> q) noexcept
{
    // std::abs
    const auto c = q._count_internal();
    return Quantity<U>(c < 0 ? -c : c);
}

template <typename U>
[[nodiscard]] Quantity<U> fabs(Quantity<U> q) noexcept
{
    return Quantity<U>(std::fabs(q._count_internal()));
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Q, typename Z>
[[nodiscard]] constexpr Absolute<Q, Z> abs(Absolute<Q, Z> q) noexcept
{
    // std::abs
    const auto c = q._count_internal();
    return Absolute<Q, Z>(c < 0 ? -c : c);
}

template <typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> fabs(Absolute<Q, Z> q) noexcept
{
    return Absolute<Q, Z>(std::fabs(q._count_internal()));
}

//==================================================================================================
// min/max
//==================================================================================================

template <typename U>
[[nodiscard]] constexpr Quantity<U> min(Quantity<U> x, Quantity<U> y) noexcept
{
    // std::min
    const auto cx = x._count_internal();
    const auto cy = y._count_internal();
    return Quantity<U>(cx < cy ? cx : cy);
}

template <typename U>
[[nodiscard]] constexpr Quantity<U> max(Quantity<U> x, Quantity<U> y) noexcept
{
    // std::max
    const auto cx = x._count_internal();
    const auto cy = y._count_internal();
    return Quantity<U>(cx < cy ? cy : cx);
}

template <typename U>
[[nodiscard]] Quantity<U> fmin(Quantity<U> x, Quantity<U> y) noexcept
{
    return Quantity<U>(std::fmin(x._count_internal(), y._count_internal()));
}

template <typename U>
[[nodiscard]] Quantity<U> fmax(Quantity<U> x, Quantity<U> y) noexcept
{
    return Quantity<U>(std::fmax(x._count_internal(), y._count_internal()));
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Q, typename Z>
[[nodiscard]] constexpr Absolute<Q, Z> min(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    // std::min
    const auto cx = x._count_internal();
    const auto cy = y._count_internal();
    return Absolute<Q, Z>(cx < cy ? cx : cy);
}

template <typename Q, typename Z>
[[nodiscard]] constexpr Absolute<Q, Z> max(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    // std::max
    const auto cx = x._count_internal();
    const auto cy = y._count_internal();
    return Absolute<Q, Z>(cx < cy ? cy : cx);
}

template <typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> fmin(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    return Absolute<Q, Z>(std::fmin(x._count_internal(), y._count_internal()));
}

template <typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> fmax(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    return Absolute<Q, Z>(std::fmax(x._count_internal(), y._count_internal()));
}

//==================================================================================================
// pow
//==================================================================================================

// q^N
template <int64_t N, typename U>
[[nodiscard]] auto pow(Quantity<U> q)
{
    using conversion    = typename U::conversion;   // Conversion<Ratio, exp>
    using kind          = typename U::kind;         // Kind<dimension, Tag>
    using dimension     = typename U::dimension;    // Ratio

    static_assert(N >= 1,
        "only positive exponents allowed");
    static_assert(std::is_same_v<typename kind::tag, kinds::Untagged>,
        "operation not supported - instead of 'pow<N>(q)' you must use 'pow<N>(q.untagged())'");
    //
    // TODO:
    // Complex-tags?
    //

    using conversion_pow = impl::RationalPower<typename conversion::ratio, N>;
    using dimension_pow = impl::RationalPower<dimension, N>;

    using C = typename Conversion<typename conversion_pow::ratio, conversion::exp * N>::type;
    using D = typename dimension_pow::ratio;

    using K = typename Kind<D, kinds::Untagged>::type;

    if constexpr (N == 1)
    {
        return Quantity<typename Unit<C, K>::type>(q._count_internal());
    }
    else if constexpr (N == 2)
    {
        return Quantity<typename Unit<C, K>::type>(q._count_internal() * q._count_internal());
    }
    else
    {
        return Quantity<typename Unit<C, K>::type>(std::pow(q._count_internal(), static_cast<Scalar>(N)));
    }
}

// q^2
template <typename U>
[[nodiscard]] auto square(Quantity<U> q)
{
    return uom::pow<2>(q);
}

// q^3
template <typename U>
[[nodiscard]] auto cube(Quantity<U> q)
{
    return uom::pow<3>(q);
}

//==================================================================================================
// nth_root
//==================================================================================================

// q^(1/N)
template <int64_t N, typename U>
[[nodiscard]] auto nth_root(Quantity<U> q)
{
    using conversion    = typename U::conversion;   // Conversion<Ratio, exp>
    using kind          = typename U::kind;         // Kind<dimension, Tag>
    using dimension     = typename U::dimension;    // Ratio

    static_assert(N >= 1,
        "only positive exponents allowed");
    static_assert(std::is_same_v<typename kind::tag, kinds::Untagged>,
        "operation not supported - instead of 'nth_root<N>(q)' you must use 'nth_root<N>(q.untagged())'");
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

    using C = typename Conversion<typename conversion_root::ratio, conversion::exp / N>::type;
    using D = typename dimension_root::ratio;

    using K = typename Kind<D, kinds::Untagged>::type;

    if constexpr (N == 1)
    {
        return Quantity<typename Unit<C, K>::type>(q._count_internal());
    }
    else if constexpr (N == 2)
    {
        return Quantity<typename Unit<C, K>::type>(std::sqrt(q._count_internal()));
    }
    else if constexpr (N == 3)
    {
        return Quantity<typename Unit<C, K>::type>(std::cbrt(q._count_internal()));
    }
    else
    {
        return Quantity<typename Unit<C, K>::type>(std::pow(q._count_internal(), 1 / static_cast<Scalar>(N)));
    }
}

// q^(1/2)
template <typename U>
[[nodiscard]] auto sqrt(Quantity<U> q)
{
    return uom::nth_root<2>(q);
}

// q^(1/3)
template <typename U>
[[nodiscard]] auto cbrt(Quantity<U> q)
{
    return uom::nth_root<3>(q);
}

//==================================================================================================
// hypot
//==================================================================================================

template <typename U>
[[nodiscard]] Quantity<U> hypot(Quantity<U> x, Quantity<U> y) noexcept
{
    return Quantity<U>(std::hypot(x._count_internal(), y._count_internal()));
}

template <typename U>
[[nodiscard]] Quantity<U> hypot(Quantity<U> x, Quantity<U> y, Quantity<U> z) noexcept
{
    return Quantity<U>(std::hypot(x._count_internal(), y._count_internal(), z._count_internal()));
}

//==================================================================================================
// exp
//==================================================================================================

[[nodiscard]] inline Dimensionless exp(Dimensionless q) noexcept
{
    return Dimensionless(std::exp(q._count_internal()));
}

[[nodiscard]] inline Dimensionless exp2(Dimensionless q) noexcept
{
    return Dimensionless(std::exp2(q._count_internal()));
}

[[nodiscard]] inline Dimensionless expm1(Dimensionless q) noexcept
{
    return Dimensionless(std::expm1(q._count_internal()));
}

[[nodiscard]] inline Dimensionless log(Dimensionless q) noexcept
{
    return Dimensionless(std::log(q._count_internal()));
}

[[nodiscard]] inline Dimensionless log10(Dimensionless q) noexcept
{
    return Dimensionless(std::log10(q._count_internal()));
}

[[nodiscard]] inline Dimensionless log1p(Dimensionless q) noexcept
{
    return Dimensionless(std::log1p(q._count_internal()));
}

[[nodiscard]] inline Dimensionless log2(Dimensionless q) noexcept
{
    return Dimensionless(std::log2(q._count_internal()));
}

//==================================================================================================
// Trigonometric
//==================================================================================================

#if 0
[[nodiscard]] inline Dimensionless sin(Dimensionless x) noexcept
{
    return Dimensionless(std::sin(x._count_internal()));
}

[[nodiscard]] inline Dimensionless cos(Dimensionless x) noexcept
{
    return Dimensionless(std::cos(x._count_internal()));
}

[[nodiscard]] inline Dimensionless tan(Dimensionless x) noexcept
{
    return Dimensionless(std::tan(x._count_internal()));
}

[[nodiscard]] inline Dimensionless asin(Dimensionless x) noexcept
{
    return Dimensionless(std::asin(x._count_internal()));
}

[[nodiscard]] inline Dimensionless acos(Dimensionless x) noexcept
{
    return Dimensionless(std::acos(x._count_internal()));
}

[[nodiscard]] inline Dimensionless atan(Dimensionless x) noexcept
{
    return Dimensionless(std::atan(x._count_internal()));
}
#else
[[nodiscard]] inline Dimensionless sin(Radians x) noexcept
{
    return Dimensionless(std::sin(x._count_internal()));
}

[[nodiscard]] inline Dimensionless cos(Radians x) noexcept
{
    return Dimensionless(std::cos(x._count_internal()));
}

[[nodiscard]] inline Dimensionless tan(Radians x) noexcept
{
    return Dimensionless(std::tan(x._count_internal()));
}

[[nodiscard]] inline Radians asin(Dimensionless x) noexcept
{
    return Radians(std::asin(x._count_internal()));
}

[[nodiscard]] inline Radians acos(Dimensionless x) noexcept
{
    return Radians(std::acos(x._count_internal()));
}

[[nodiscard]] inline Radians atan(Dimensionless x) noexcept
{
    return Radians(std::atan(x._count_internal()));
}

template <typename U>
[[nodiscard]] Radians atan2(Quantity<U> y, Quantity<U> x) noexcept
{
    return Radians(std::atan2(y._count_internal(), x._count_internal()));
}
#endif

[[nodiscard]] inline Dimensionless sinh(Dimensionless x) noexcept
{
//  return (exp(x) - exp(-x)) / 2;
    return Dimensionless(std::sinh(x._count_internal()));
}

[[nodiscard]] inline Dimensionless cosh(Dimensionless x) noexcept
{
//  return (exp(x) + exp(-x)) / 2;
    return Dimensionless(std::cosh(x._count_internal()));
}

[[nodiscard]] inline Dimensionless tanh(Dimensionless x) noexcept
{
//  return sinh(x) / cosh(x);
    return Dimensionless(std::tanh(x._count_internal()));
}

[[nodiscard]] inline Dimensionless asinh(Dimensionless x) noexcept
{
    return Dimensionless(std::asinh(x._count_internal()));
}

[[nodiscard]] inline Dimensionless acosh(Dimensionless x) noexcept
{
    return Dimensionless(std::acosh(x._count_internal()));
}

[[nodiscard]] inline Dimensionless atanh(Dimensionless x) noexcept
{
    return Dimensionless(std::atanh(x._count_internal()));
}

//==================================================================================================
//
//==================================================================================================

namespace impl
{
    template <typename Scale, typename U, typename Fn>
    auto ToInt(Quantity<U> x, Fn fn) noexcept
    {
        static_assert(IsRatio<Scale> || IsQuantity<Scale>,
            "operation not supported");

        if constexpr (IsRatio<Scale>)
        {
            using T = ScaledQuantity<Conversion<Scale>, Quantity<U>>;
            return ( T( fn( x.as(T{})._count_internal() ) ) ).as(Quantity<U>{});
        }
        else
        {
            return Scale( fn( x.as(Scale{})._count_internal() ) );
        }
    }

    template <typename Scale, typename Q, typename Z, typename Fn>
    auto ToInt(Absolute<Q, Z> x, Fn fn) noexcept
    {
        static_assert(IsRatio<Scale>,
            "operation not supported");

        return Absolute<Q, Z>( ToInt<Scale>( Q(x._count_internal()), fn )._count_internal() );
    }
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Scale = Ratio<1>, typename U>
[[nodiscard]] auto ceil(Quantity<U> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const Scalar t) { return std::ceil(t); } );
}

template <typename Scale = Ratio<1>, typename U>
[[nodiscard]] auto floor(Quantity<U> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const Scalar t) { return std::floor(t); } );
}

template <typename Scale = Ratio<1>, typename U>
[[nodiscard]] auto trunc(Quantity<U> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const Scalar t) { return std::trunc(t); } );
}

template <typename Scale = Ratio<1>, typename U>
[[nodiscard]] auto round(Quantity<U> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const Scalar t) { return std::round(t); } );
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Scale = Ratio<1>, typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> ceil(Absolute<Q, Z> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const Scalar t) { return std::ceil(t); } );
}

template <typename Scale = Ratio<1>, typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> floor(Absolute<Q, Z> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const Scalar t) { return std::floor(t); } );
}

template <typename Scale = Ratio<1>, typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> trunc(Absolute<Q, Z> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const Scalar t) { return std::trunc(t); } );
}

template <typename Scale = Ratio<1>, typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> round(Absolute<Q, Z> x) noexcept
{
    return impl::ToInt<Scale>(x, [](const Scalar t) { return std::round(t); } );
}

//==================================================================================================
//
//==================================================================================================

template <typename U1, typename U2, typename U3, std::enable_if_t<std::is_same_v<uom::impl::MulUnits<U1, U2>, U3>, int> = 0>
[[nodiscard]] auto fma(Quantity<U1> x, Quantity<U2> y, Quantity<U3> z)
{
    // The function signature is more restrictive than (x * y + z).
    // This is to avoid any implicit conversion taking place, which would defeat the purpose
    // the fma instruction.

    return Quantity<U3>(std::fma(x._count_internal(), y._count_internal(), z._count_internal()));
}

//==================================================================================================
// midpoint
//==================================================================================================

template <typename U>
[[nodiscard]] Quantity<U> midpoint(Quantity<U> x, Quantity<U> y) noexcept
{
    return Quantity<U>(x._count_internal() / 2 + y._count_internal() / 2);
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> midpoint(Absolute<Q, Z> x, Absolute<Q, Z> y) noexcept
{
    return Absolute<Q, Z>(x._count_internal() / 2 + y._count_internal() / 2);
}

//==================================================================================================
// lerp
//==================================================================================================

template <typename U>
[[nodiscard]] Quantity<U> lerp(Quantity<U> x, Quantity<U> y, Scalar t) noexcept
{
    return Quantity<U>((1 - t) * x._count_internal() + t * y._count_internal());
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Q, typename Z>
[[nodiscard]] Absolute<Q, Z> lerp(Absolute<Q, Z> x, Absolute<Q, Z> y, Scalar t) noexcept
{
    return Absolute<Q, Z>((1 - t) * x._count_internal() + t * y._count_internal());
}

//==================================================================================================
//
//==================================================================================================

namespace impl
{
//  inline constexpr Scalar kPi    = static_cast<Scalar>(3.14159265358979323846264338328);
    inline constexpr Scalar kTwoPi = static_cast<Scalar>(6.28318530717958647692528676656);

    inline Scalar Clamp(Scalar x, Scalar lo, Scalar hi)
    {
        if (x < lo)
            return lo;
        if (x > hi)
            return hi;

        return x;
    }

    inline Scalar NormalizePositive(Scalar x, Scalar range) noexcept
    {
        // UNITS_ASSERT(std::isfinite(x));

        Scalar w = x;
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

    inline Scalar NormalizeSymmetric(Scalar x, Scalar range) noexcept
    {
        // UNITS_ASSERT(std::isfinite(x));

        const Scalar half = range / 2;

        Scalar w = x;
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
[[nodiscard]] inline Radians normalize_positive(Radians x) noexcept
{
    return Radians(impl::NormalizePositive(x._count_internal(), impl::kTwoPi));
}

// Reduce x to [-pi, pi)
[[nodiscard]] inline Radians normalize_symmetric(Radians x) noexcept
{
    return Radians(impl::NormalizeSymmetric(x._count_internal(), impl::kTwoPi));
}

// Reduce x to [0, 360)
[[nodiscard]] inline Degrees normalize_positive(Degrees x) noexcept
{
    return Degrees(impl::NormalizePositive(x._count_internal(), 360));
}

// Reduce x to [-180, 180)
[[nodiscard]] inline Degrees normalize_symmetric(Degrees x) noexcept
{
    return Degrees(impl::NormalizeSymmetric(x._count_internal(), 360));
}

// Reduce x to [0, 400)
[[nodiscard]] inline Gons normalize_positive(Gons x) noexcept
{
    return Gons(impl::NormalizePositive(x._count_internal(), 400));
}

// Reduce x to [-200, 200)
[[nodiscard]] inline Gons normalize_symmetric(Gons x) noexcept
{
    return Gons(impl::NormalizeSymmetric(x._count_internal(), 400));
}

// Reduce x to [0, 1)
[[nodiscard]] inline Revolutions normalize_positive(Revolutions x) noexcept
{
    return Revolutions(impl::NormalizePositive(x._count_internal(), 1));
}

// Reduce x to [-1/2, 1/2)
[[nodiscard]] inline Revolutions normalize_symmetric(Revolutions x) noexcept
{
    return Revolutions(impl::NormalizeSymmetric(x._count_internal(), 1));
}

} // namespace uom
